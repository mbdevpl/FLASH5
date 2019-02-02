!!****if* source/Driver/DriverMain/Driver_verifyInitDt
!!
!! NAME
!!  Driver_verifyInitDt
!!
!! SYNOPSIS
!!  Driver_verifyInitDt()
!!
!! DESCRIPTION
!!
!! The initial timestep "dt" is a runtime parameter for the simulations.
!! This routine makes sure that users haven't inadvertently provided
!! the initial value for dt that violates the Courant-Friedrichs-Levy
!! condition.
!!
!! ARGUMENTS
!!
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

#include "Flash.h"
#include "constants.h" 
#undef FIXEDBLOCKSIZE
subroutine Driver_verifyInitDt()

  use Driver_data, ONLY : dr_restart, dr_dt, dr_dtInit, dr_dtOld, dr_globalMe,&
       dr_dtSTS, dr_dtNew, dr_meshComm,                                       &
       dr_globalComm,dr_dtDiffuse, dr_dtAdvect,dr_dtHeatExch,dr_useSTS,       &
       dr_tstepSlowStartFactor
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use Hydro_interface, ONLY : Hydro_computeDt, Hydro_consolidateCFL
  use Diffuse_interface, ONLY: Diffuse_computeDt
  use flash_iterator, ONLY : flash_iterator_t
  use flash_tile, ONLY : flash_tile_t

  implicit none       

  include "Flash_mpi.h"
  

  real,dimension(3) :: dtCheck  ,dtCFL

  integer    :: dtMinLoc(5)
  integer :: i, ierr

  integer :: coordSize
  logical :: gcell = .true.
  real, dimension(MDIM) :: del
#ifdef INDEXREORDER
  integer,parameter::IX=1,IY=2,IZ=3
#else
  integer,parameter::IX=2,IY=3,IZ=4
#endif  

  real, allocatable,dimension(:)::&
       xCoord,dx,uxgrid,yCoord,dy,uygrid,zCoord,dz,uzgrid,&
       xLeft,xRight,yLeft,yRight,zLeft,zRight

  !arrays which hold the starting and ending indicies of a block
  integer,dimension(2,MDIM)::lim

  !!coordinate infomration to be passed into physics  
  real, pointer :: solnData(:,:,:,:)
  integer :: isize,jsize,ksize
  logical :: runVerifyInitDt = .false.
  real :: extraHydroInfo
  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc
  integer:: ib

  nullify(solnData)

!!$  dr_dtSTS = 0.0     !First use is in a max(dr_dtSTS,...), see Driver_evolveFlash. - KW
!!$  dr_dtNew = 0.0     !First use is in a max(dr_dtSTS,...), see Driver_evolveFlash. - KW

  !! Need to run this routine when the super time stepping algorithm is used.
  if (dr_useSTS) then
     runVerifyInitDt = .true.
  else
     if (.not. dr_restart) then
        runVerifyInitDt = .true.
     endif
  endif

  if (.not. dr_restart) then
     ! compute the CFL timestep for the simulation and compare it to the
     ! user specified initial timestep.  Scream loudly if there is a problem.
     
     !initialize values 
     dtCheck = huge(dtCheck)
     dtMinLoc(:) = 0
     
#ifdef CFL_VAR
     call Hydro_consolidateCFL()
#endif
        
     !There is some overhead in calling Hydro_computeDt.  Although it is a
     !pain to get the coordinates and solution data before calling the 
     !routine, this is just initialization.  Getting the coordinates inside
     !Hydro_computeDt would be much more costly during the run

     call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
      
        call tileDesc%getDataPtr(solnData, CENTER)

        ! DEV: Do we really need to get the data for the GC as well?
        lim(:, :) = tileDesc%limits(:, :)
        lim(LOW,  1:NDIM) = lim(LOW,  1:NDIM) - NGUARD
        lim(HIGH, 1:NDIM) = lim(HIGH, 1:NDIM) + NGUARD
        associate(lo => lim(LOW,  :), &
                  hi => lim(HIGH, :))
           allocate(xCoord(lo(IAXIS):hi(IAXIS)))
           allocate(dx(    lo(IAXIS):hi(IAXIS)))
           allocate(uxgrid(lo(IAXIS):hi(IAXIS)))
           allocate(yCoord(lo(JAXIS):hi(JAXIS)))
           allocate(dy(    lo(JAXIS):hi(JAXIS)))
           allocate(uygrid(lo(JAXIS):hi(JAXIS)))
           allocate(zCoord(lo(KAXIS):hi(KAXIS)))
           allocate(dz(    lo(KAXIS):hi(KAXIS)))
           allocate(uzgrid(lo(KAXIS):hi(KAXIS)))
           allocate(xLeft( lo(IAXIS):hi(IAXIS)))
           allocate(xRight(lo(IAXIS):hi(IAXIS)))
           allocate(yLeft( lo(JAXIS):hi(JAXIS)))
           allocate(yRight(lo(JAXIS):hi(JAXIS)))
           allocate(zLeft( lo(KAXIS):hi(KAXIS)))
           allocate(zRight(lo(KAXIS):hi(KAXIS)))
        end associate 

        call tileDesc%coordinates(IAXIS, CENTER,     TILE_AND_HALO, xCoord)
        call tileDesc%coordinates(IAXIS, LEFT_EDGE,  TILE_AND_HALO, xLeft)
        call tileDesc%coordinates(IAXIS, RIGHT_EDGE, TILE_AND_HALO, xRight)
        
        if (NDIM > 1) then
           call tileDesc%coordinates(JAXIS, CENTER,     TILE_AND_HALO, yCoord)
           call tileDesc%coordinates(JAXIS, LEFT_EDGE,  TILE_AND_HALO, yLeft)
           call tileDesc%coordinates(JAXIS, RIGHT_EDGE, TILE_AND_HALO, yRight)

           if (NDIM > 2) then
              call tileDesc%coordinates(KAXIS, CENTER,     TILE_AND_HALO, zCoord)
              call tileDesc%coordinates(KAXIS, LEFT_EDGE,  TILE_AND_HALO, zLeft)
              call tileDesc%coordinates(KAXIS, RIGHT_EDGE, TILE_AND_HALO, zRight)
           endif
        endif

        call tileDesc%deltas(del)
        dx(:) = del(1)
        dy(:) = del(2)
        dz(:) = del(3)
        
        uxgrid(:) = 0
        uygrid(:) = 0
        uzgrid(:) = 0
        
        call Hydro_computeDt (tileDesc, &
             xCoord, dx, uxgrid, &
             yCoord, dy, uygrid, &
             zCoord, dz, uzgrid, &
             tileDesc%limits, lim,  &
             solnData,      &
             dtCheck(1), dtMinLoc, &
             extraInfo=extraHydroInfo)
!!$     call Diffuse_computeDt ( blockList(i), &
!!$          xCoord, xLeft,xRight, dx, uxgrid, &
!!$          yCoord, yLeft,yRight, dy, uygrid, &
!!$          zCoord, zLeft,zRight, dz, uzgrid, &
!!$          lim,limGC,  &
!!$          solnData,      &
!!$          dtCheck(2), dtMinLoc )
        
        call tileDesc%releaseDataPtr(solnData, CENTER)
 
        deallocate(xCoord)
        deallocate(dx)
        deallocate(uxgrid)
        deallocate(yCoord)
        deallocate(dy)
        deallocate(uygrid)
        deallocate(zCoord)
        deallocate(dz)
        deallocate(uzgrid)
        deallocate(xLeft)
        deallocate(xRight)
        deallocate(yLeft)
        deallocate(yRight)
        deallocate(zLeft)
        deallocate(zRight)
        
        call itor%next()
     end do
     call Grid_releaseTileIterator(itor)
     
     ! find the minimum across all processors, store it in dtCFL on MasterPE
     call MPI_AllReduce(dtCheck(1), dtCFL(1), 3, FLASH_REAL, MPI_MIN, &
          dr_globalComm, ierr)
     
     
     !! Initialize advection and diffusion time steps
     dr_dtAdvect  = dtCFL(1)
     dr_dtDiffuse = dtCFL(2)
     dr_dtHeatExch= dtCFL(3)
     
     if (.not. dr_useSTS) then
        dtCFL(1) = minval(dtCFL)
     endif
     
     if (dr_dtInit > dr_tstepSlowStartFactor*dtCFL(1)) then
        
        if (dr_globalMe .EQ. MASTER_PE) then
           print *, '***********************************************************'
           print *, ' Warning: The initial timestep is too large.'
           print *, '   initial timestep = ', dr_dtInit
           print *, '   CFL timestep     = ', dtCFL(1)
           print *, ' Resetting dtinit to dr_tstepSlowStartFactor*dtcfl.'
           print *, '***********************************************************'
           print *, ' '
        endif
        
        dr_dt = dr_tstepSlowStartFactor*dtCFL(1)
        
     else
        
        dr_dt = dr_dtInit
        
     endif
     
     dr_dtOld = dr_dt
     !print *, dr_dt, "dr_dt initial final"
     !!     call Cosmology_computeDt(dtCheck(1))
     
  else
     
     dr_dtAdvect = dr_dt
     dr_dtDiffuse = dr_dt
     
  endif
  return
end subroutine Driver_verifyInitDt








