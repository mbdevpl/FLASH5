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
  use Grid_interface, ONLY :  Grid_getCellCoords, Grid_getDeltas, &
       Grid_getSingleCellCoords, Grid_getMaxRefinement, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getLeafIterator, Grid_releaseLeafIterator
  use Hydro_interface, ONLY : Hydro_computeDt, Hydro_consolidateCFL
  use Diffuse_interface, ONLY: Diffuse_computeDt
  use leaf_iterator, ONLY : leaf_iterator_t
  use block_metadata, ONLY : block_metadata_t

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
  integer,dimension(2,MDIM)::lim,limGC

  !!coordinate infomration to be passed into physics  
  real, pointer :: solnData(:,:,:,:)
  integer :: isize,jsize,ksize
  logical :: runVerifyInitDt = .false.
  real :: extraHydroInfo
  type(leaf_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc
  integer:: ib, level, maxLev

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

  call Grid_getMaxRefinement(maxLev,mode=1)
  
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
        
     do level=1,maxLev
        call Grid_getLeafIterator(itor, level=level)
        do while(itor%is_valid())
           call itor%blkMetaData(blockDesc)
          
           lim=blockDesc%Limits
           limGC=blockDesc%LimitsGC
           call Grid_getBlkPtr(blockDesc, solnData)

           isize = limGC(HIGH,IAXIS)-limGC(LOW,IAXIS)+1
           jsize = limGC(HIGH,JAXIS)-limGC(LOW,JAXIS)+1
           ksize = limGC(HIGH,KAXIS)-limGC(LOW,KAXIS)+1
           
           allocate(xCoord(limGC(LOW,IAXIS):limGC(HIGH,IAXIS)))
           allocate(dx(limGC(LOW,IAXIS):limGC(HIGH,IAXIS)))
           allocate(uxgrid(limGC(LOW,IAXIS):limGC(HIGH,IAXIS)))
           allocate(yCoord(limGC(LOW,JAXIS):limGC(HIGH,JAXIS)))
           allocate(dy(limGC(LOW,JAXIS):limGC(HIGH,JAXIS)))
           allocate(uygrid(limGC(LOW,JAXIS):limGC(HIGH,JAXIS)))
           allocate(zCoord(limGC(LOW,KAXIS):limGC(HIGH,KAXIS)))
           allocate(dz(limGC(LOW,KAXIS):limGC(HIGH,KAXIS)))
           allocate(uzgrid(limGC(LOW,KAXIS):limGC(HIGH,KAXIS)))
           allocate(xLeft(limGC(LOW,IAXIS):limGC(HIGH,IAXIS)))
           allocate(xRight(limGC(LOW,IAXIS):limGC(HIGH,IAXIS)))
           allocate(yLeft(limGC(LOW,JAXIS):limGC(HIGH,JAXIS)))
           allocate(yRight(limGC(LOW,JAXIS):limGC(HIGH,JAXIS)))
           allocate(zLeft(limGC(LOW,KAXIS):limGC(HIGH,KAXIS)))
           allocate(zRight(limGC(LOW,KAXIS):limGC(HIGH,KAXIS)))
           
           
           coordSize = isize
           call Grid_getCellCoords(IAXIS,blockDesc,CENTER,gcell,xCoord,coordSize)
           call Grid_getCellCoords(IAXIS,blockDesc,LEFT_EDGE,gcell,xLeft,isize)
           call Grid_getCellCoords(IAXIS,blockDesc,RIGHT_EDGE,gcell,xRight,isize)
           
           if (NDIM > 1) then
              coordSize = jsize
              call Grid_getCellCoords(JAXIS,blockDesc,CENTER,gcell,yCoord,coordSize)
              call Grid_getCellCoords(JAXIS,blockDesc,LEFT_EDGE,gcell,yLeft,jsize)
              call Grid_getCellCoords(JAXIS,blockDesc,RIGHT_EDGE,gcell,yRight,jsize)
              
              if (NDIM > 2) then
                 coordSize = ksize
                 call Grid_getCellCoords(KAXIS,blockDesc,CENTER,gcell,zCoord,coordSize)
                 call Grid_getCellCoords(KAXIS,blockDesc,LEFT_EDGE,gcell,zLeft,ksize)
                 call Grid_getCellCoords(KAXIS,blockDesc,RIGHT_EDGE,gcell,zRight,ksize)
              endif
           endif
           
           
           
           call Grid_getDeltas(level, del)
           dx(:) = del(1)
           dy(:) = del(2)
           dz(:) = del(3)
           
           uxgrid(:) = 0
           uygrid(:) = 0
           uzgrid(:) = 0
           
           call Hydro_computeDt (blockDesc, &
                xCoord, dx, uxgrid, &
                yCoord, dy, uygrid, &
                zCoord, dz, uzgrid, &
                lim,limGC,  &
                solnData,      &
                dtCheck(1), dtMinLoc, &
                extraInfo=extraHydroInfo)
!!$        call Diffuse_computeDt ( blockList(i), &
!!$             xCoord, xLeft,xRight, dx, uxgrid, &
!!$             yCoord, yLeft,yRight, dy, uygrid, &
!!$             zCoord, zLeft,zRight, dz, uzgrid, &
!!$             lim,limGC,  &
!!$             solnData,      &
!!$             dtCheck(2), dtMinLoc )
           
           call Grid_releaseBlkPtr(blockDesc, solnData)
           nullify(solnData)
 
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
        call Grid_releaseLeafIterator(itor)
     end do
     
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








