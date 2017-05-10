!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_Unsplit/Hydro_computeDt
!!
!! NAME
!!
!!  Hydro_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Hydro_computeDt(integer(IN)    :: blockID,
!!                  real(IN)       :: x(GRID_ILO_GC:GRID_IHI_GC),
!!                  real(IN)       :: dx(GRID_ILO_GC:GRID_IHI_GC),
!!                  real(IN)       :: uxgrid(GRID_ILO_GC:GRID_IHI_GC),
!!                  real(IN)       :: y(GRID_JLO_GC:GRID_JHI_GC),
!!                  real(IN)       :: dy(GRID_JLO_GC:GRID_JHI_GC),
!!                  real(IN)       :: uygrid(GRID_JLO_GC:GRID_JHI_GC),
!!                  real(IN)       :: z(GRID_KLO_GC:GRID_KHI_GC),
!!                  real(IN)       :: dz(GRID_KLO_GC:GRID_KHI_GC),
!!                  real(IN)       :: uzgrid(GRID_KLO_GC:GRID_KHI_GC),
!!                  integer(IN)    :: blkLimits(2,MDIM)
!!                  integer(IN)    :: blkLimitsGC(2,MDIM)
!!                  real,pointer   :: U(:,:,:,:),
!!                  real(INOUT)    :: dtCheck,
!!                  integer(INOUT) :: dtMinLoc(5),
!!                  real(INOUT), optional :: extraInfo)
!!  
!!
!! DESCRIPTION
!!
!!  This routine computes the timestep limiter for the Unsplit Hydro solver.
!!  The Courant-Fredrichs-Lewy criterion is used.  The sound
!!  speed is computed and together with the velocities, is used to constrain
!!  the timestep such that no information can propagate more than one zone
!!  per timestep.
!!  Note that this routine only accounts for computing advection time step in hyperbolic
!!  system of equations.
!!
!! ARGUMENTS
!!
!!  blockID       -  local block ID
!!  x, y, z       -  three, directional coordinates
!!  dx,dy,dz      -  distances in each {*=x, y z} directions
!!  uxgrid        -  velocity of grid expansion in x directions
!!  uygrid        -  velocity of grid expansion in y directions
!!  uzgrid        -  velocity of grid expansion in z directions
!!  blkLimits     -  the indices for the interior endpoints of the block
!!  blkLimitsGC   -  the indices for endpoints including the guardcells
!!  U             -  the physical, solution data from grid
!!  dtCheck       -  variable to hold timestep constraint
!!  dtMinLoc(5)   -  array to hold location of cell responsible for minimum dt:
!!                   dtMinLoc(1) = i index
!!                   dtMinLoc(2) = j index
!!                   dtMinLoc(3) = k index
!!                   dtMinLoc(4) = 1 if found to meet criterion, 0 otherwise
!!                   dtMinLoc(5) = hy_meshMe
!!  extraInfo     -  Driver_computeDt can provide extra info to the caller
!!                   using this argument.
!!
!!***

!!REORDER(4): U

Subroutine Hydro_computeDt( &
     x, dx, uxgrid, &
     y, dy, uygrid, &
     z, dz, uzgrid, &
     blkLimits, blkLimitsGC, &
     U,  dtCheck, dtMinLoc,  &
     extraInfo)


#include "Flash.h"
#include "constants.h"

  use Hydro_data,       ONLY : hy_geometry, hy_cfl, hy_dref, hy_eref, &
                               hy_pref, hy_vref, hy_meshMe,           &
                               hy_useHydro, hy_updateHydroFluxes,     &
                               hy_useVaryingCFL
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

  !! Arguments type declaration ------------------------------------------
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif

  real, pointer         :: U(:,:,:,:)
  real,   intent(INOUT) :: dtCheck
  integer,intent(INOUT) :: dtMinLoc(5)
  real, OPTIONAL,intent(INOUT) :: extraInfo
  !! ----------------------------------------------------------------------

  integer :: i, j, k, temploc(5)
  real    :: sndspd2, delxinv, delyinv, delzinv, dt_temp, dt_ltemp


  if ((.not. hy_useHydro) .or. (.not. hy_updateHydroFluxes)) return

  dt_temp    = 0.
  temploc(:) = 0

  !! Conversion for unitSystem if needed ---------------------------------------
  U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
  U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
  U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
  U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)/hy_vref
  !! ---------------------------------------------------------------------------


  if (NDIM == 1) then
     !--------------------------------------------------------------------------!
     ! 1-dimensional                                                            !
     !--------------------------------------------------------------------------!
     delxinv = 1.0/dx(blkLimits(LOW,IAXIS))

     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

        sndspd2   = U(GAMC_VAR,i,1,1)*U(PRES_VAR,i,1,1)/U(DENS_VAR,i,1,1) 
        dt_ltemp  = (abs(U(VELX_VAR,i,1,1)-uxgrid(i))+sqrt(sndspd2))*delxinv

        if (dt_ltemp > dt_temp) then
           dt_temp    = dt_ltemp
           temploc(1) = i
           temploc(2) = 1
           temploc(3) = 1
           temploc(4) = 1
           temploc(5) = hy_meshMe
        endif

     enddo

  elseif (NDIM == 2) then
     !--------------------------------------------------------------------------!
     ! 2-dimensional                                                            !
     !--------------------------------------------------------------------------!
     if (hy_geometry == CARTESIAN .OR. hy_geometry == CYLINDRICAL) then 
        ! the 'y' coordinate is not angular in Cartesian and cylindrical coords

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))

        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              sndspd2  = U(GAMC_VAR,i,j,1)*U(PRES_VAR,i,j,1)/U(DENS_VAR,i,j,1)

              dt_ltemp = max((abs(U(VELX_VAR,i,j,1)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                   (abs(U(VELY_VAR,i,j,1)-uygrid(j))+sqrt(sndspd2))*delyinv)

              if (dt_ltemp > dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = 1
                 temploc(5) = hy_meshMe
              endif

           enddo
        enddo
               
     else ! Angular coordinates in 2D: Spherical or Polar

        ! y is angular
        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))

        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              delyinv = 1.0/(x(i)*dy(j))

              sndspd2  = U(GAMC_VAR,i,j,1)*U(PRES_VAR,i,j,1)/U(DENS_VAR,i,j,1)
              dt_ltemp = max((abs(U(VELX_VAR,i,j,1)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                   (abs(U(VELY_VAR,i,j,1)-uygrid(j))+sqrt(sndspd2))*delyinv)

              if (dt_ltemp > dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = 1
                 temploc(5) = hy_meshMe
              endif
           enddo
        enddo

     endif

  elseif (NDIM == 3) then
     !--------------------------------------------------------------------------!
     ! 3-dimensional                                                            !
     !--------------------------------------------------------------------------!

     if (hy_geometry == CARTESIAN) then

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))
        delzinv = 1.0/dz(blkLimits(LOW,KAXIS))

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 sndspd2  = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 dt_ltemp = max((abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                      (abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(sndspd2))*delyinv,&
                      (abs(U(VELZ_VAR,i,j,k)-uzgrid(k))+sqrt(sndspd2))*delzinv)

                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = 1
                    temploc(5) = hy_meshMe
                 endif

              enddo
           enddo
        enddo

     elseif (hy_geometry == CYLINDRICAL) then

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 delzinv = 1.0/(x(i)*dz(k)) ! z is phi

                 sndspd2  = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 dt_ltemp = max((abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                      (abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(sndspd2))*delyinv,&
                      (abs(U(VELZ_VAR,i,j,k)-uzgrid(k))+sqrt(sndspd2))*delzinv)

                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = 1
                    temploc(5) = hy_meshMe
                 endif

              enddo
           enddo
        enddo
        
     elseif (hy_geometry == SPHERICAL) then

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 delyinv = 1.0/(x(i)*dy(j)) ! y is theta
                 delzinv = 1.0/(x(i)*sin(y(j))*dz(k)) ! z is phi

                 sndspd2  = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 dt_ltemp = max((abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                      (abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(sndspd2))*delyinv,&
                      (abs(U(VELZ_VAR,i,j,k)-uzgrid(k))+sqrt(sndspd2))*delzinv)

                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = 1
                    temploc(5) = hy_meshMe
                 endif

              enddo
           enddo
        enddo

     else ! Polar in 3D (that's a no no)
        call Driver_abortFlash("[Hydro_computeDt] ERROR: Polar geometry not supported in 3D")
     endif

  endif


  dt_temp = hy_cfl / dt_temp
  if (dt_temp < dtCheck) then
     dtCheck = dt_temp
     dtMinLoc = temploc
  endif

  !! For the purpose of having screen output for CFL
  if (present(extraInfo)) then
     if (hy_useVaryingCFL) then
        extraInfo = hy_cfl
     else
        extraInfo = 0.
     end if
  end if

  !! Conversion for unitSystem if needed ---------------------------------------
  U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)*hy_dref
  U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)*hy_eref
  U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)*hy_pref
  U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)*hy_vref
  !! ---------------------------------------------------------------------------
  if(dtCheck <= 0.0) call Driver_abortFlash("[Hydro]: Computed dt is not positive! Aborting!")
  return

End Subroutine Hydro_computeDt


