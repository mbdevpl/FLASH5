!!****if* source/Simulation/SimulationMain/MacLaurin/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Maclaurin spheroid problem.
!!
!!  References:  Maclaurin, C. 1742, 
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!  blockID     --   The number of the block to initialize
!!
!!***
!!REORDER(4):solnData
subroutine Simulation_initBlock(solnData,block)

  use Simulation_data, ONLY: sim_density, sim_gamma, sim_Omega2, sim_Pconst, &
       sim_nsubinv, sim_nsubzones, sim_xctr, sim_yctr, sim_zctr, &
       sim_initGeometry, sim_geom2DAxisymmetric, sim_geom3DCartesian, &
       sim_a3inv, sim_a1inv, &
       sim_smallRho, sim_smallP, sim_smallE
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getDeltas
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"

  real,                   pointer    :: solnData(:,:,:,:)
  type(block_metadata_t), intent(in) :: block

  integer, dimension(LOW:HIGH,MDIM) :: blkLimitsGC
  integer, dimension(MDIM) :: startingPos
  real, dimension(MDIM) :: del
  integer  :: sizeX, sizeY, sizeZ
  logical  :: gcell=.true.
  integer  :: i, j, k, ii, jj, kk
  real     :: xdist, ydist, zdist, dist2, rxy, rxyz2, rinv
  real     :: xx, yy, zz, dxx, dyy, dzz, vxfac, vyfac, vzfac
  real     :: sum_rho, sum_p, sum_vx, sum_vy, sum_vz, pres, vel

  real, dimension(:), allocatable :: xl, yl, zl
  real  :: vx, vy, vz, p, rho, e, ek, ei, gam


  ! Get the coordinate information for the current block

  ! get the coordinate information for the current block
  blkLimitsGC = block%limitsGC
  
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
  allocate(xl(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
  allocate(yl(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
  allocate(zl(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, block, LEFT_EDGE, gcell, zl, sizeZ)
  endif
  if (NDIM >= 2) then 
     call Grid_getCellCoords(JAXIS, block, LEFT_EDGE, gcell, yl, sizeY)
  endif
  call Grid_getCellCoords(IAXIS, block, LEFT_EDGE, gcell, xl, sizeX)
  
  ! Set initial conditions in each zone
  call Grid_getDeltas(block%level,del)
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     dzz = del(KAXIS) * sim_nsubinv
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        dyy = del(JAXIS) * sim_nsubinv
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
           dxx = del(IAXIS) * sim_nsubinv

           sum_rho = 0.0
           sum_p   = 0.0
           sum_vx  = 0.0
           sum_vy  = 0.0
           sum_vz  = 0.0

           ! Break the zone into nsubzones^ndim sub-zones and average the
           ! results to get the values for the zone.

           do kk = 0, sim_nsubzones-1
              zz    = zl(k) + (kk + 0.5)*dzz 
              zdist = (zz - sim_zctr) * K3D

              do jj = 0, sim_nsubzones-1
                 yy    = yl(j) + (jj + 0.5)*dyy
                 ydist = (yy - sim_yctr) * K2D

                 do ii = 0, sim_nsubzones-1
                    xx    = xl(i) + (ii + 0.5)*dxx
                    xdist = xx - sim_xctr

                    select case (sim_initGeometry)

                    case (sim_geom2DAxisymmetric)

                       dist2 = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                       rxy   = xdist
                       rxyz2 = dist2
                       rinv  = 1./sqrt(xdist**2 + ydist**2)
                       vxfac = 0.0
                       vyfac = 0.0
                       vzfac = 0.0

                    case (sim_geom3DCartesian)

                       dist2 = (xdist*sim_a1inv)**2 + (ydist*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                       rxy   = sqrt(xdist**2 + ydist**2)
                       rxyz2 = (rxy*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                       rinv  = 1./sqrt(xdist**2 + ydist**2 + zdist**2)
                       vxfac = -ydist*rinv
                       vyfac = xdist*rinv
                       vzfac = 0.0

                    end select

                    if (dist2 <= 1.) then    ! inside the spheroid
                       pres    = sim_Pconst * (1.0 - rxyz2)
                       vel     = rxy * sim_Omega2
                       sum_rho = sum_rho + sim_density
                       sum_p   = sum_p   + pres
                       sum_vx  = sum_vx  + vel*vxfac
                       sum_vy  = sum_vy  + vel*vyfac
                       sum_vz  = sum_vz  + vel*vzfac
                    else                     ! outside the spheroid
                       sum_rho = sum_rho + sim_smallRho*10
                       sum_p   = sum_p   + sim_smallP*10
                    endif

                 enddo
              enddo
           enddo

           rho = max(sum_rho * sim_nsubinv**3, sim_smallRho)
           p   = max(sum_p * sim_nsubinv**3, sim_smallP)
           vx  = sum_vx*sim_nsubinv**3
           vy  = sum_vy*sim_nsubinv**3
           vz  = sum_vz*sim_nsubinv**3
           ek  = 0.5*(vx**2+vy**2+vz**2)
           gam = sim_gamma
           ei  = max(p/(rho*(gam-1.0)), sim_smallE)
           e   = ei + ek
           solnData(DENS_VAR,i,j,k) =  rho
           solnData(PRES_VAR,i,j,k) =  p
           solnData(ENER_VAR,i,j,k) =  e
           solnData(GAME_VAR,i,j,k) =  gam
           solnData(GAMC_VAR,i,j,k) =  gam
           solnData(VELX_VAR,i,j,k) =  vx
           solnData(VELY_VAR,i,j,k) =  vy
           solnData(VELZ_VAR,i,j,k) =  vz
           solnData(EINT_VAR,i,j,k) =  ei
           
        enddo
     enddo
  enddo
  deallocate(xl, yl, zl)
  return
end subroutine Simulation_initBlock

