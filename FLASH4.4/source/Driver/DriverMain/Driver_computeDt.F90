!!****if* source/Driver/DriverMain/Driver_computeDt
!!
!! NAME
!!
!!  Driver_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Driver_computeDt(integer(IN) :: nbegin,
!!                  integer(IN) :: nstep, 
!!                  real(IN)    :: simTime,
!!                  real(IN)    :: dtOld,
!!                  real(OUT)   :: dtNew)
!!
!! DESCRIPTION
!!
!!  Determine the stability-limited time step.
!!  This timestep is determined using information from the included
!!  physics modules - many different timestep limiters are polled.
!!
!!  The global driver might use a different (hopefully smaller) time
!!  step, to match a file write time (tplot or trstr) or if the
!!  simulation end time has been reached; such possibilities are
!!  not considered here.
!!
!! ARGUMENTS
!!
!!  nbegin - first step of the simulation (nbegin is only used
!!              to determine if a label header should be written to
!!              the screen)
!!  nstep - current step of the simulation
!!  simTime - current simulation time of the run
!!  dtOld - the dt from the timestep that we just finished 
!!         (it's old because we be using dtOld to calculate 
!!          and return the dt for the next timestep (dtNew)
!!  dtNew - returned value of the dt calculated for the next timestep
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!! The calls to units currently not included in the code are commented out.
!!
!!
!!
!!*** 
#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif
#include "constants.h"
#include "Flash.h"
#undef FIXEDBLOCKSIZE
subroutine Driver_computeDt(nbegin, nstep, &
                    simTime, dtOld, dtNew)

  use Driver_data, ONLY : dr_dtMin,dr_dtMax, dr_tstepChangeFactor, &
                          dr_redshift, dr_useRedshift,             &
                          dr_printTStepLoc,                        &
                          dr_dtSTS, dr_useSTS, dr_globalMe, dr_globalComm,&
                          dr_dtAdvect, dr_dtDiffuse,dr_meshComm,     &
                          dr_dtMinContinue, dr_dtMinBelowAction
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface,ONLY : Logfile_stamp
  use IO_interface,     ONLY : IO_writeCheckpoint
  use Grid_interface, ONLY : Grid_getSingleCellCoords,Grid_getMaxRefinement, &
       Grid_getTileIterator, Grid_releaseTileIterator
  use Hydro_interface, ONLY : Hydro_computeDt, Hydro_consolidateCFL
  use Heat_interface, ONLY : Heat_computeDt
  use Diffuse_interface, ONLY : Diffuse_computeDt 
  use Burn_interface, ONLY : Burn_computeDt
  use RadTrans_interface, ONLY: RadTrans_computeDt
  use Particles_interface, ONLY: Particles_computeDt

  use IncompNS_interface, ONLY : IncompNS_computeDt
  use flash_iterator, ONLY : flash_iterator_t
  use flash_tile, ONLY : flash_tile_t 

  implicit none

  include "Flash_mpi.h"

  integer, intent(IN) :: nbegin, nstep
  real,    intent(IN) :: simTime    !! current simulation time
  real, intent(IN) :: dtOld      !! last time step we used
  real, intent(OUT):: dtNew      !! the new timestep we get. to be returned.
 

  ! Local variables and functions
  integer :: i, error, iout

  !! This is arbitrarily fixed. Users that add more units that influence the time
  !! should change this.

  integer, parameter :: nUnits = 15
  real, PARAMETER :: MAX_TSTEP = huge(1.0)

  
  real    :: dtModule(2,nUnits), dtLocal(2,nUnits)
  integer :: dtMinLoc(5), lminloc(5,nUnits), ngmin, pgmin
  integer :: status(MPI_Status_Size)

  logical :: gcell = .true.
  real, DIMENSION(MDIM) :: coords

  real, dimension(nUnits) :: tstepOutput
  character (len=20), save, DIMENSION(nUnits) :: &
                         limiterName, limiterNameOutput

  !!prepatory data structures for passing coords to timestep routines
  real, dimension(MDIM) :: del
  integer, dimension(MDIM) :: index

#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_ILO_GC:GRID_IHI_GC) :: xLeft,xRight,xCenter
  real,dimension(GRID_JLO_GC:GRID_JHI_GC) :: yLeft,yRight,yCenter
  real,dimension(GRID_KLO_GC:GRID_KHI_GC) :: zLeft,zRight,zCenter
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) ::  dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) ::  dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) ::  dz, uzgrid
#else
  real, allocatable,dimension(:)::&
       dx,uxgrid,dy,uygrid,dz,uzgrid
  real, allocatable,dimension(:)::xLeft,xRight,xCenter
  real, allocatable,dimension(:)::yLeft,yRight,yCenter
  real, allocatable,dimension(:)::zLeft,zRight,zCenter

#endif

  !arrays which hold the starting and ending indicies of a block
  integer,dimension(LOW:HIGH,MDIM)::blkLimitsGC

  !!coordinate infomration to be passed into physics  
  real, pointer :: solnData(:,:,:,:)
  integer :: isize,jsize,ksize

  logical :: printTStepLoc
  integer :: itempLimit = 0
  integer, parameter :: HYDRO=1,BURN=2,GRAV=3,HEAT=4,COOL=5,TEMP=6,&
                        PART=7,DIFF=8,COSMO=9,STIR=10,HEATXCHG=11, &
                        RADTRANS=12,STS=13,INS=14,SOLID=15
#ifdef INDEXREORDER
  integer,parameter::IX=1,IY=2,IZ=3
#else
  integer,parameter::IX=2,IY=3,IZ=4
#endif  
  logical, save :: firstCall = .TRUE.
  logical :: limitersChanged
  logical :: printToScrn
  logical :: endRun
  real :: extraHydroInfo
  character (len=10) :: cflNumber
  character (len=10) :: extraHydroInfoStr
  character (len=10*nUnits) :: tailStr
  character (len=24) :: tailFmt
  real :: extraHydroInfoMin
  real :: extraHydroInfoApp
  real :: dtNewComputed
  type(flash_iterator_t) :: itor
  type(flash_tile_t) :: tileDesc
  integer:: ib, level, maxLev
  real :: err

  nullify(solnData)

  ! Initializing extraHydroInfo to zero:
  extraHydroInfo = 0.
  extraHydroInfoMin = 1.e10 !temporary large fake CFL for comparison
  extraHydroInfoApp = 1.e10 !temporary large fake CFL


  data limiterName(HYDRO) /'dt_hydro'/
  data limiterName(HEAT) /'dt_Heat'/
  data limiterName(PART) /'dt_Part '/
  data limiterName(BURN) /'dt_Burn '/
  data limiterName(COOL) /'dt_Cool '/
  data limiterName(TEMP) /'dt_Temp '/
  data limiterName(DIFF) /'dt_Diff '/
  data limiterName(COSMO) /'dt_Cosm'/
  data limiterName(STIR) /'dt_Stir'/
  data limiterName(HEATXCHG) /'dt_HeatXchg'/
  data limiterName(RADTRANS) /'dt_RadTrans'/
  data limiterName(STS)  /'dt_STS'/
  data limiterName(INS)  /'dt_INS'/
  data limiterName(SOLID) /'dt_Solid'/
  data limiterNameOutput /nUnits * ' '/
  data cflNumber  /'CFL'/


  !            Find the local minimum timestep among the included physics
  !            modules for locally stored blocks.
  
  !            Initialize all timestep variables.




  printTStepLoc = dr_printTStepLoc
  endRun = .FALSE.
  
  dtMinLoc(:) = 0
  lminloc(:,:) = 0
  lminloc(NDIM+1:MDIM,:) = 1

  do i = 1, nUnits
     dtLocal(1,i) = MAX_TSTEP
     dtLocal(2,i) = real(dr_globalMe)
  enddo

  !            Loop over all local leaf-node blocks
  
  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.

  call Hydro_consolidateCFL()
  do level=1,maxLev
     call Grid_getTileIterator(itor, LEAF, level=level, tiling=.FALSE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)

        ! Match the size of the FIXEDBLOCKSIZE case
        ! This disallows tiling.
        blkLimitsGC = tileDesc%blkLimitsGC
        call tileDesc%getDataPtr(solnData, CENTER)

#ifndef FIXEDBLOCKSIZE
        associate(lo => blkLimitsGC(LOW,  :), &
                  hi => blkLimitsGC(HIGH, :))
           allocate(xLeft  (lo(IAXIS):hi(IAXIS)))
           allocate(xRight (lo(IAXIS):hi(IAXIS)))
           allocate(xCenter(lo(IAXIS):hi(IAXIS)))
           allocate(dx     (lo(IAXIS):hi(IAXIS)))
           allocate(uxgrid (lo(IAXIS):hi(IAXIS)))
           allocate(yLeft  (lo(JAXIS):hi(JAXIS)))
           allocate(yRight (lo(JAXIS):hi(JAXIS)))
           allocate(yCenter(lo(JAXIS):hi(JAXIS)))
           allocate(dy     (lo(JAXIS):hi(JAXIS)))
           allocate(uygrid (lo(JAXIS):hi(JAXIS)))
           allocate(zLeft  (lo(KAXIS):hi(KAXIS)))
           allocate(zRight (lo(KAXIS):hi(KAXIS)))
           allocate(zCenter(lo(KAXIS):hi(KAXIS)))
           allocate(dz     (lo(KAXIS):hi(KAXIS)))
           allocate(uzgrid (lo(KAXIS):hi(KAXIS)))
        end associate
#endif
#ifdef DEBUG_DRIVER
        print*,'before calling get coordinates'
#endif
        call tileDesc%coordinates(IAXIS, CENTER,    TILE_AND_HALO, xCenter)
        call tileDesc%coordinates(IAXIS, LEFT_EDGE, TILE_AND_HALO, xLeft)
        call tileDesc%coordinates(IAXIS, RIGHT_EDGE,TILE_AND_HALO, xRight)
        
#ifdef DEBUG_DRIVER
        print*,'before calling get coordinates'
#endif
        if (NDIM > 1) then
           call tileDesc%coordinates(JAXIS, CENTER,    TILE_AND_HALO, yCenter)
           call tileDesc%coordinates(JAXIS, LEFT_EDGE, TILE_AND_HALO, yLeft)
           call tileDesc%coordinates(JAXIS, RIGHT_EDGE,TILE_AND_HALO, yRight)
           
           if (NDIM > 2) then
#ifdef DEBUG_DRIVER
              print*,'before calling get coordinates'
#endif
              call tileDesc%coordinates(KAXIS, CENTER,    TILE_AND_HALO, zCenter)
              call tileDesc%coordinates(KAXIS, LEFT_EDGE, TILE_AND_HALO, zLeft)
              call tileDesc%coordinates(KAXIS, RIGHT_EDGE,TILE_AND_HALO, zRight)
           endif
        endif
        
        call tileDesc%deltas(del)
        dx(:) = del(1)
        dy(:) = del(2)
        dz(:) = del(3)
        
        uxgrid(:) = 0
        uygrid(:) = 0
        uzgrid(:) = 0
        
        ! hydro
#ifdef DEBUG_DRIVER
        print*,'going to call Hydro timestep'
#endif
        !extraHydroInfo = 0.
        call Hydro_computeDt (tileDesc, &
             xCenter, dx, uxgrid, &
             yCenter, dy, uygrid, &
             zCenter, dz, uzgrid, &
             tileDesc%limits, blkLimitsGC,  &
             solnData,      &
             dtLocal(1,HYDRO), lminloc(:,HYDRO), &
             extraInfo=extraHydroInfo )
        
        !! Extra CFL information
        if (extraHydroInfo .ne. 0.) then
           if (extraHydroInfo <= extraHydroInfoMin) then
              extraHydroInfoMin = extraHydroInfo
           endif
           if (lminloc(4,HYDRO) == level) then
              extraHydroInfoApp = extraHydroInfo
           endif
        else !if extraHydroInfo == 0.
           extraHydroInfoMin = 0.0
           extraHydroInfoApp = 0.0
        endif
        
        
#ifdef DEBUG_DRIVER
        print*,'returned from hydro timestep'
#endif
        
#ifndef FIXEDBLOCKSIZE
        deallocate(xCenter)
        deallocate(xLeft)
        deallocate(xRight)
        deallocate(dx)
        deallocate(uxgrid)
        deallocate(yCenter)
        deallocate(yLeft)
        deallocate(yRight)
        deallocate(dy)
        deallocate(uygrid)
        deallocate(zCenter)
        deallocate(zLeft)
        deallocate(zRight)
        deallocate(dz)
        deallocate(uzgrid)
#endif
#ifdef DEBUG_DRIVER
        print*,'release blockpointer'
#endif
        
        call tileDesc%releaseDataPtr(solnData, CENTER)
        nullify(solnData)

        call itor%next()
     enddo
     call Grid_releaseTileIterator(itor)
  end do
!!$     end do
     
!!$  !! Choose the smallest CFL for screen output - provisional, may change below
  extraHydroInfo = 0.
  call MPI_AllReduce (extraHydroInfoMin, extraHydroInfo, 1, & 
       FLASH_REAL, MPI_MIN, dr_globalComm, error)


!!$  ! IncompNS:
!!$  call IncompNS_computeDt(dtLocal(1,INS),lminloc(:,INS))


  ! DEV: we disabled temperature timestep limiter for now.  
  ! The old temperature was not updated with the refinement, 
  ! so dT/T was precomputed and the number of blocks may not be
  ! the same as there are now.
!!$  if (itempLimit == 1) then
!!$     do blockID = 1,  MAXBLOCKS
!!$       call Driver_computeDtTemp(dr_globalMe, dtOld, dtLocal(1,6), &
!!$             lminloc(1,6), blockID)
!!$     enddo
!!$  endif
  


  !            Find the minimum timestep across all processors and all
  !            modules.

  call MPI_AllReduce (dtLocal(1,1), dtModule(1,1), nUnits, & 
       MPI_2Double_Precision, MPI_MinLoc, dr_globalComm, error)

  dtNew    = huge(dtNew)    ! dt will hold the minimum timestep
  ngmin = 1                 ! ngmin will hold the winning module #
  pgmin = MASTER_PE         ! pgmin will hold the winning PE #
  
!!$  do i = 1, nUnits-1
!!$     if (dtModule(1,i) < dtNew) then
!!$        dtNew = dtModule(1,i)
!!$        pgmin = dtModule(2,i)
!!$        ngmin = i
!!$     endif
!!$  enddo

  do i = 1, nUnits
     if ((i .ne. STS) .and. (i .ne. DIFF)) then
        if (dtModule(1,i) < dtNew) then
           dtNew = dtModule(1,i)
           pgmin = dtModule(2,i)
           ngmin = i
        endif
     endif
  enddo

  ! Save it to hydro's advection time scale
  ! Note: This dr_dtAdvect is the minimum timestep from all physics units,
  !       e.g., Hydro, Stir, Burn, Heat, Cool, Particle, and Cosmology,
  !       except from DIFF and STS.
  dr_dtAdvect = dtNew

  ! Do it one more time
  if (dtModule(1,DIFF) < dtNew) then
     dtNew = dtModule(1,DIFF)
     pgmin = dtModule(2,DIFF)
     ngmin = DIFF
  endif

  ! Save it to hydro's diffusion time scale
  dr_dtDiffuse = dtModule(1,DIFF)
      
  ! have the processor that is determining the timestep limit broadcast the
  ! proc number, block number, and i,j,k of the zone that set the timestep
  ! to all processors

  dtMinLoc(:) = lminloc(:,ngmin)

  call MPI_Bcast(dtMinLoc(1), 5, MPI_INTEGER, pgmin, dr_globalComm, error)

  if (extraHydroInfo .NE. 0.0 .AND. ngmin == HYDRO) then
  ! If the unit that determines the timestep is HYDRO, have the
  ! processor that is determining the timestep limit broadcast the
  ! local cfl of the zone that set the timestep (which should now be
  ! in extraHydroInfoApp) to all processors
     extraHydroInfo = extraHydroInfoApp
     call MPI_Bcast(extraHydroInfo, 1, FLASH_REAL, pgmin, dr_globalComm, error)
  end if

  if (dtNew < dr_dtMinContinue) then
     endRun = .TRUE.
     dtNewComputed = dtNew
     if (dr_globalMe == MASTER_PE) then
        print*,'        About to exit, computed time step is too small:',dtNew
     end if
  end if
     

  ! limit the timestep to increase by at most a factor of dr_tstepChangeFactor

  dtNew = min( dtNew, dtOld*dr_tstepChangeFactor )
  
  if (nstep .GE. 50) then   !! This is where Cellular starts to fail
!     print *, 'nstep = ',nstep
  endif

  ! Use dr_dtmin and dr_dtmax to limit the timestep.  If this makes the code
  ! unstable, it's not our fault.
  dtNew = min( max( dtNew, dr_dtMin ), dr_dtMax )


  if (printTStepLoc) then
         
     ! convert the dtMinLoc array into a physical position (x,y,z) where the
     ! timestep is being set.  dtMinLoc(5) is the processor number, dtMinLoc(4)
     ! is the blockID on that proc.
     coords(:) = 0.0

     if (dr_globalMe == dtMinLoc(5)) then

        if (dtMinLoc(4) > 0) then
           index(:)=dtMinLoc(1:MDIM)
           call Grid_getSingleCellCoords(index,dtMinLoc(4),CENTER,coords=coords)
        else
           coords(:) = 999.0
        end if

        ! send this to the master processor
        if (dr_globalMe /= MASTER_PE) then

           call MPI_Send (coords(1), 3, FLASH_REAL, MASTER_PE, & 
                0, dr_globalComm, error)
           
        endif
        
     elseif (dr_globalMe == MASTER_PE) then
        
        call MPI_Recv (coords(1), 3, FLASH_REAL, dtMinLoc(5), 0, & 
             dr_globalComm, status, error)            
        
     endif
     
  endif
  
  ! Print out the time and next timestep.
  
  ! only print out the timestep from the limiters that are active
  iout = 0
  limitersChanged = .FALSE.
  do i = 1, nUnits
     if (dtModule(1,i) /= MAX_TSTEP) then
        iout = iout + 1
        if (.NOT.(firstCall.OR.limitersChanged)) then
           if (limiterNameOutput(iout) .NE. limiterName(i)) limitersChanged = .TRUE.
        end if
        tstepOutput(iout) = dtModule(1,i)
        limiterNameOutput(iout) = limiterName(i)
     endif
  enddo

!!$print*,iout,nUnits;pause


  printToScrn = .true.
  if (printToScrn .AND. dr_globalMe == MASTER_PE) then


     if (extraHydroInfo .eq. 0.) then

        if (printTStepLoc) then
        
           if (nstep == nbegin .OR. limitersChanged) then

              if (.not. dr_useRedshift) then
                 write (*,803) 'n', 't', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout)
              else
                 write (*,804) 'n', 't', 'z', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout)
              endif
              
           endif
        
           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,801) nstep, simTime, dtNew, coords(1), coords(2), &
                      coords(3), (tstepOutput(i),i=1,iout)
              else
                 write(*,801) nstep, simTime, max(dtNew,dr_dtSTS), coords(1), coords(2), &
                      coords(3), (tstepOutput(i),i=1,iout)
              endif

           else
              if (.not. dr_useSTS) then
                 write(*,802) nstep, simTime, dr_redshift, dtNew, coords(1), &
                      coords(2), coords(3), (tstepOutput(i),i=1,iout)
              else
                 write(*,802) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), coords(1), &
                      coords(2), coords(3), (tstepOutput(i),i=1,iout)
              endif
           endif
           
        else
        
           if (nstep .eq. nbegin .OR. limitersChanged) then
           
              if (.not. dr_useRedshift) then
                 write (*,903) 'n', 't', 'dt', (limiterNameOutput(i),i=1,iout)
              else
                 write (*,904) 'n', 't', 'z', 'dt', (limiterNameOutput(i),i=1,iout)
              endif
           
           endif
        
           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,901) nstep, simTime, dtNew, (tstepOutput(i),i=1,iout)
              else
                 write(*,901) nstep, simTime, max(dtNew,dr_dtSTS), (tstepOutput(i),i=1,iout)
              endif
           else
              if (.not. dr_useSTS) then
                 write(*,902) nstep, simTime, dr_redshift, dtNew, &
                      (tstepOutput(i),i=1,iout)
              else
                 write(*,902) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), &
                      (tstepOutput(i),i=1,iout)
              endif

           endif

        endif

     else ! elseif (extraHydroInfo .ne. 0.) then

        if (printTStepLoc) then
           if (nstep == nbegin .OR. limitersChanged) then

              if (.not. dr_useRedshift) then
                 write (*,803) 'n', 't', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout),cflNumber
              else
                 write (*,804) 'n', 't', 'z', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout),cflNumber
              endif
              
           endif
        
           if (extraHydroInfo .GE. 0.0001 .AND. extraHydroInfo .LT. 9.9995) then
              write(extraHydroInfoStr,807) extraHydroInfo
           else
              write(extraHydroInfoStr,808) extraHydroInfo
           end if
!!$           write (tailFmt,'(A1,I2,A)') "(",iout,"(:,1X,ES9.3),1x,A10)"
!!$           write (tailStr,tailFmt) (tstepOutput(i),i=1,iout), adjustl(extraHydroInfoStr)

           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,805) nstep, simTime, dtNew, coords(1), coords(2), &
                      coords(3), trim(tailStr)
              else
                 write(*,805) nstep, simTime, max(dtNew,dr_dtSTS), coords(1), coords(2), &
                      coords(3), trim(tailStr)
              endif

           else
              if (.not. dr_useSTS) then
                 write(*,806) nstep, simTime, dr_redshift, dtNew, coords(1), &
                      coords(2), coords(3), trim(tailStr)
              else
                 write(*,806) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), coords(1), &
                      coords(2), coords(3), trim(tailStr)
              endif
           endif
           
        else
        
           if (nstep .eq. nbegin .OR. limitersChanged) then
           
              if (.not. dr_useRedshift) then
                 write (*,903) 'n', 't', 'dt', (limiterNameOutput(i),i=1,iout),cflNumber
              else
                 write (*,904) 'n', 't', 'z', 'dt', (limiterNameOutput(i),i=1,iout),cflNumber
              endif
           
           endif
        
           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,901) nstep, simTime, dtNew, (tstepOutput(i),i=1,iout), extraHydroInfo
              else
                 write(*,901) nstep, simTime, max(dtNew,dr_dtSTS), (tstepOutput(i),i=1,iout), extraHydroInfo
              endif
           else
              if (.not. dr_useSTS) then
                 write(*,902) nstep, simTime, dr_redshift, dtNew, &
                      (tstepOutput(i),i=1,iout), extraHydroInfo
              else
                 write(*,902) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), &
                      (tstepOutput(i),i=1,iout), extraHydroInfo
              endif

           endif

        endif

     endif


  endif ! end of (printToScrn .and. dr_globalMe == MASTER_PE)
      

801 format (1X, I7, 1x, ES10.4, 1x, ES10.4, 2x, '(', 1P,G10.3, ', ', &
            G 10.3, ', ', G10.3, ')', ' | ', 11(:, 1X, ES9.3))
802 format (1X, I7, 1x, ES10.4, 1x, F8.3, 1x, ES10.4, 2x, '(', ES9.3, ', ', &
            ES 9.3, ', ', ES9.3, ')', ' | ', 11(:, 1X, ES9.3))
803 format (1X, A7, 1x, A10, 1x, A10, 2x, '(', A10, ', ', A10, ', ', A10, ')', &
            ' | ', 11(:, 1X, A9),1x,A10)
804 format (1X, A7, 1x, A10, 1x, A7, 1x, A10, 2x, '(', A9, ', ', A9, ', ', &
         A9, ')', ' | ', 11(:, 1X, A9),1x,A10)
805 format (1X, I7, 1x, ES10.4, 1x, ES10.4, 2x, '(', 1P,G10.3, ', ', &
            G 10.3, ', ', G10.3, ')', ' | ', A)
806 format (1X, I7, 1x, ES10.4, 1x, F8.3, 1x, ES10.4, 2x, '(', ES9.3, ', ', &
            ES 9.3, ', ', ES9.3, ')', ' | ', A)
807 format (F10.7)
808 format (ES10.3)

901 format (1X, I7, 1X, ES10.4, 1x, ES10.4, ' | ', 11(:, 1X, ES11.5),1x,ES10.4)
902 format (1X, I7, 1X, ES10.4, 1x, F8.3, 1x, ES10.4, ' | ', 11(:, 1X, ES11.5),1x,ES10.4)
903 format (1X, A7, 1x, A10   , 1x, A10,    ' | ', 11(:, 1X, A11),1x,A10)
904 format (1X, A7, 1x, A10   , 1x, A7, 1x, A10,    ' | ', 11(:, 1X, A11),1x,A10)

  if (endRun) then
     call Logfile_stamp(dtNewComputed, '[Driver_computeDt] Computed dtNew')
     call Logfile_stamp(dtNew        , '[Driver_computeDt] Next dtNew would be')
     if (dr_dtMinBelowAction==1) then
        call Logfile_stamp( 'Writing additional checkpoint because of dr_dtMinBelowAction' , '[Driver_computeDt]')
        call IO_writeCheckpoint()
     end if
     call Logfile_stamp( 'Exiting simulation because dr_dtNew < dr_dtMinContinue' , '[Driver_computeDt]')
     call Driver_abortFlash('Computed new time step smaller than dr_dtMinContinue!')
  end if

  firstCall = .FALSE.
  return
end subroutine Driver_computeDt


!!$     call Burn_computeDt ( blockID,  &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,               &
!!$                           dtLocal(1,BURN), lminloc(:,BURN) )
           
!!$     call Gravity_computeDt ( blockID, dr_globalMe, &
!!$                           xCenter,xLeft,xRight, dx, uxgrid, &
!!$                           yCenter,yLeft,yRight, dy, uygrid, &
!!$                           zCenter,zLeft,zRight, dz, uzgrid, &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,      &
!!$                           dtLocal(1,GRAV), lminloc(:,GRAV) )
     
     
!!$     call Heat_computeDt ( blockID,  &
!!$                           xCenter, dx, uxgrid, &
!!$                           yCenter, dy, uygrid, &
!!$                           zCenter, dz, uzgrid, &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,      &
!!$                           dtLocal(1,HEAT), lminloc(:,HEAT) )
!!$
!!$     call RadTrans_computeDt(blockID, blkLimits,blkLimitsGC, &
!!$          solnData, dtLocal(1,RADTRANS), lminloc(:,RADTRANS) )
!!$
!!$     call Particles_computeDt &
!!$          ( blockID, dtLocal(1,PART), lminloc(:,PART))
!!$
!!$
!!$     call Diffuse_computeDt ( blockID,  &
!!$                           xCenter,xLeft,xRight, dx, uxgrid, &
!!$                           yCenter,yLeft,yRight, dy, uygrid, &
!!$                           zCenter,zLeft,zRight, dz, uzgrid, &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,      &
!!$                           dtLocal(1,DIFF), lminloc(:,DIFF) )
!!$     

!!$      !! Super time step
!!$      if (dr_useSTS) then
!!$         dtLocal(1,STS) = dr_dtSTS
!!$      endif

