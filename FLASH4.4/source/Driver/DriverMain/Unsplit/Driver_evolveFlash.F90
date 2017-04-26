!!****if* source/Driver/DriverMain/Unsplit/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  call Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This routine implements a directionally unsplit integrator for time 
!! advancement. A single step in this driver accounts for all 
!! multidimensional fluxes and updates conserved variables for the next
!! time step. The routine is a default driver for the Unsplit Staggered
!! Mesh (USM) MHD and Unsplit Hydro solvers.
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_",
!! like dr_globalMe, dr_dt, or dr_beginStep, and are stored in Fortran
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_".
!!
!!
!!***

!#define DEBUG_DRIVER
#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif


subroutine Driver_evolveFlash()

  use Driver_data,         ONLY : dr_globalMe, dr_globalNumProcs, dr_nbegin, &
                                  dr_nend, dr_dt,                        &
                                  dr_tmax, dr_simTime, dr_redshift,      &
                                  dr_nstep, dr_dtOld, dr_dtNew,          &
                                  dr_simGeneration,                      &
                                  dr_restart,                            &
                                  dr_redshiftOld, dr_useRedshift,        &
                                  dr_redshiftfinal,                      &
                                  dr_useSTS, dr_nuSTS, dr_nstepTotalSTS, &
                                  dr_dtSTS,dr_dt_subSTS,                 &
                                  dr_dtAdvect, dr_dtDiffuse,             &
                                  dr_useSTSforDiffusion,                 &
                                  dr_tstepChangeFactor,                  &
                                  dr_allowDtSTSDominate,dr_meshComm
  use Driver_interface,    ONLY : Driver_sourceTerms, Driver_computeDt, &
                                  Driver_superTimeStep, &
                                  Driver_logMemoryUsage, &
                                  Driver_driftUnk, &
                                  Driver_diagnostics
  use Logfile_interface,   ONLY : Logfile_stamp, Logfile_close
  use Timers_interface,    ONLY : Timers_start, Timers_stop, &
                                  Timers_getSummary
  use Diffuse_interface,   ONLY : Diffuse
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface,      ONLY : Grid_updateRefinement,&
                                  Grid_fillGuardCells,&
                                  Grid_releaseBlkPtr
  use Hydro_interface,     ONLY : Hydro, &
                                  Hydro_gravPotIsAlreadyUpdated
  use Gravity_interface,   ONLY : Gravity_potentialListOfBlocks
  use IO_interface,        ONLY : IO_output,IO_outputFinal
  use RadTrans_interface,  ONLY : RadTrans
  use Eos_interface,       ONLY : Eos_logDiagnostics
  use Simulation_interface, ONLY: Simulation_adjustEvolution
  use Profiler_interface, ONLY : Profiler_start, Profiler_stop
  use famrex_multivab_module, ONLY: famrex_multivab, famrex_multivab_build, &
                                    famrex_mviter, famrex_mviter_build,&
                                    famrex_mviter_destroy,famrex_multivab_destroy
  use famrex_box_module,      ONLY: famrex_box


  implicit none

#include "constants.h"
#include "Flash.h"

  integer   :: localNumBlocks

  integer,save :: sweepDummy = SWEEP_ALL

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(4,2) :: strBuff
  character(len=15) :: numToStr

  logical :: gridChanged
  logical :: endRunPl !Should we end our run on this iteration, based on conditions detected by the IO unit?
  logical :: endRun !Should we end our run on this iteration, based on conditions detected by the IO unit?
  logical :: endRunWallClock !Should we end our run on this iteration, based on wall clock time?

  ! for super-time-stepping
  integer :: nstepSTS
  real    :: dt_diffuse_temp
  real    :: dtNewTemp
  logical :: useSTS_local
  integer :: nstepTotalSTS_local 

  integer, parameter :: driftUnk_flags = DRIFT_NO_PARENTS
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif
  integer:: ib, blockID
  integer, dimension(LOW:HIGH,MDIM) :: tileLimits,blkLimitsGC
  integer :: blockCount
  integer,dimension(MAXBLOCKS)::blks
  real,pointer,dimension(:,:,:,:) :: Uout
  real,dimension(MDIM) :: del

  type(famrex_multivab),target :: phi
  type(famrex_mviter) :: mvi
  type(famrex_box) :: bx, tbx


  endRunPl = .false.
  endRun = .false.

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Profiler_start("FLASH_evolution")
  call Timers_start("evolution")

  do dr_nstep = dr_nBegin, dr_nend
     
     useSTS_local = dr_useSTS

     if (dr_globalMe == MASTER_PE) then
        
        write (numToStr(1:), '(I10)') dr_nstep
        write (strBuff(1,1), "(A)") "n"
        write (strBuff(1,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_simTime
        write (strBuff(2,1), "(A)") "t"
        write (strBuff(2,2), "(A)") trim(adjustl(numToStr))
        
        if (.not. dr_useSTS) then
           write (numToStr(1:), "(1PE12.6)") dr_dt
        else
           write (numToStr(1:), "(1PE12.6)") max(dr_dt,dr_dtSTS)
        endif
        write (strBuff(3,1), "(A)") "dt"
        write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))
        
        call Logfile_stamp( strBuff(1:3,:), 3, 2, "step")
        
        
     end if
     
!!     call Simulation_adjustEvolution(blockCount, blockList, dr_nstep, dr_dt, dr_simTime)
     
     ! 1. Cosmology-Friedmann Eqn.
     call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
     
     dr_simTime = dr_simTime + dr_dt
     dr_simGeneration = 0
     
     ! 2. Hydro/MHD
#ifdef DEBUG_DRIVER
     print*,'going into Hydro/MHD'  ! DEBUG
     print*,'going into hydro myPE=',dr_globalMe
#endif
     !!ChageForAMRex -- Here is where we put in the iterator and extract the relevant metadata
     !!ChageForAMRex -- from the iterator and then use the case statement to transfer control to the
     !!ChageForAMRex -- right implementation.
     
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        !!        call Logfile_stampVarMask(hy_gcMask, .FALSE., '[hy_hllUnsplit]', 'gcNeed')
     end if
#endif
     
     !! Guardcell filling routine
!!$     call Grid_fillGuardCells(CENTER,ALLDIR,&
!!$          maskSize=hy_gcMaskSize, mask=hy_gcMask,makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskLogged)
     
     call Grid_fillGuardCells(CENTER,ALLDIR)
     call Timers_start("Hydro")

  call famrex_multivab_build(phi, LEAF, CENTER, dr_meshComm, NUNK_VARS)
  call famrex_mviter_build(mvi, phi, tiling=.true.) !tiling is currently ignored...
  do while(mvi%next())
       bx = mvi%tilebox()

       Uout => phi%dataptr(mvi)
       tileLimits(LOW, :) = bx%lo
       tileLimits(HIGH,:) = bx%hi


!!$     call Grid_getListOfBlocks(LEAF,blks,blockCount)
!!$     do ib=1,blockCount
!!$        blockID=blks(ib)
!!$        call Grid_getBlkIndexLimits(blockID,tileLimits,blkLimitsGC,CENTER)
!!$        call Grid_getBlkPtr(blockID,Uout,CENTER)

       blockID = mvi%localIndex() !Are we cheating here?
     
       call Grid_getDeltas(blockID,del)
       
       call Hydro(del,tileLimits,Uout,dr_simTime, dr_dt, dr_dtOld,  sweepDummy)
       call Grid_releaseBlkPtr(blockID,Uout,CENTER)
    end do
    call Timers_stop("Hydro")
    call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
#ifdef DEBUG_DRIVER
     print*, 'return from Hydro/MHD timestep'  ! DEBUG
     print*,'returning from hydro myPE=',dr_globalMe
#endif
     
     
!!$     ! 8. Diagnostics
!!$     call Timers_start("diagnostics")
!!$     call Driver_diagnostics(blockCount, blockList, dr_dt)
!!$     call Timers_stop("diagnostics")
!!$#ifdef DEBUG_DRIVER
!!$     print*, 'return from Diagnostics '  ! DEBUG
!!$#endif
     
     call famrex_multivab_destroy(phi)
     call famrex_mviter_destroy(mvi)
     !! save for old dt
     dr_dtOld = dr_dt
     
     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------
     
     !output a plotfile before the grid changes
     call Timers_start("IO_output")
     
     if (.not. useSTS_local) then
        call IO_output(dr_simTime, &
             dr_dt, dr_nstep+1, dr_nbegin, endRunPl, PLOTFILE_AND_PARTICLEFILE)
     else
        call IO_output(dr_simTime, &
             dr_dtSTS, dr_nstep+1, dr_nbegin, endRunPl, PLOTFILE_AND_PARTICLEFILE)
     endif
     call Timers_stop("IO_output")
     
     
     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")
     
     if (gridChanged) dr_simGeneration = dr_simGeneration + 1
     
     ! backup needed old
     if (.not. useSTS_local) dr_dtOld = dr_dt
     
#ifdef DEBUG_DRIVER
     print*, 'going into Driver_computeDt '  ! DEBUG
#endif
     ! calculate new
     call Timers_start("Driver_computeDt")
     call Driver_computeDt(dr_nbegin,  dr_nstep,      &
          dr_simTime, dr_dtOld, dr_dtNew)
     call Timers_stop("Driver_computeDt")
#ifdef DEBUG_DRIVER
     print*, 'return from Driver_computeDt '  ! DEBUG
#endif
     
     ! store new
     if (.not. useSTS_local) dr_dt = dr_dtNew
     
     call Timers_start("IO_output")
     if (.not. useSTS_local) then
        call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun,&
             CHECKPOINT_FILE_ONLY)
     else
        call IO_output(dr_simTime,dr_dtSTS,dr_nstep+1,dr_nbegin,endRun,&
             CHECKPOINT_FILE_ONLY)
     endif
     call Timers_stop("IO_output")
     endRun = (endRunPl .OR. endRun)
     
     call Eos_logDiagnostics(.FALSE.)
     
     
     !!*****************************************************************************
     !!  Evolution Loop -- check termination conditions
     !!*****************************************************************************
     
     !Exit if a .dump_restart or .kill was found during the last step
     if(endRun) exit
     
     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift  
     !!        (also called redshiftFinal) 
     !!  (iii) the wall clock time is greater than the maximum 
     !!        (wall_clock_time_max)
     
     
     if (dr_simTime >= dr_tmax) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max SimTime"
        endif
        exit
     end if
     
     call dr_wallClockLimitExceeded(endRunWallClock)
     if (endRunWallClock) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        endif
        exit
     end if
     
  enddo
  !The value of dr_nstep after the loop is (dr_nend + 1) if the loop iterated for
  !the maximum number of times.  However, we need to retain the value that
  !dr_nstep had during the last loop iteration, otherwise the number for nstep
  !that will be stored in a final checkpoint file will be wrong.
  dr_nstep = min(dr_nstep,dr_nend)
  
  !!******************************************************************************
  !! End of Evolution Loop
  !!******************************************************************************
  
  call Timers_stop("evolution")
  call Profiler_stop("FLASH_evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  !if a file termination, this may already be done.
  if(.NOT.endRun) call IO_outputFinal()
  call Timers_getSummary( max(0,dr_nstep-dr_nbegin+1))
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()
  
  return
  
end subroutine Driver_evolveFlash
