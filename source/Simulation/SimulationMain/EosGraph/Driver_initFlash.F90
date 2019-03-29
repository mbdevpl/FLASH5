!!****if* source/Simulation/SimulationMain/EosGraph/Driver_initFlash
!!
!! NAME
!!  Driver_initFlash
!!
!! SYNOPSIS
!!
!!   Driver_initFlash()
!!
!! DESCRIPTION
!!
!!  Performs Flash initializations, which includes:
!!
!!  Call all 'init' routines in units.  Order does matter,
!!  particularly when restarting from a checkpoint file.
!!
!!  For the most part, Driver_initFlash calls other units' init
!!  routines directly, like call IO_init or call Grid_init.  This
!!  routine also makes calls to other Driver initialization routines
!!  like Driver_initMaterialProperties or Driver_initSourceTerms.
!!  These routines then call the unit-specific initialization 
!!  routines.  This level of abstraction was added to simplify
!!  the initialization calls.
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
!! HISTORY 
!!
!!  2008-03-14 KW   Moved material properties initialization up.
!!
!!***


subroutine Driver_initFlash()
  
  use Driver_data, ONLY: dr_globalComm, dr_globalMe, dr_globalNumProcs, dr_nbegin, &
       dr_initialSimTime, dr_elapsedWCTime, &
       dr_initialWCTime, dr_restart, dr_dtInit, dr_redshift,dr_particlesInitialized


  use Driver_interface, ONLY : Driver_init, &
    Driver_initMaterialProperties, Driver_initSourceTerms, &
    Driver_verifyInitDt, Driver_abortFlash, Driver_setupParallelEnv, &
    Driver_initNumericalTools
  use RuntimeParameters_interface, ONLY : RuntimeParameters_init, RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_init
  use PhysicalConstants_interface, ONLY : PhysicalConstants_init
  use Timers_interface, ONLY : Timers_init, Timers_start, Timers_stop

  use Grid_interface, ONLY : Grid_init, Grid_initDomain, &
    Grid_getListOfBlocks
  use Multispecies_interface, ONLY : Multispecies_init
  use Particles_interface, ONLY : Particles_init,  Particles_initData, &
       Particles_initForces
 
  use Eos_interface, ONLY : Eos_init, Eos_logDiagnostics
  use Hydro_interface, ONLY : Hydro_init
  use Simulation_interface, ONLY : Simulation_init
  use IO_interface, ONLY :IO_init, IO_outputInitial
  implicit none       
  
#include "constants.h"
#include "Flash.h"

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  logical :: updateRefine

  dr_elapsedWCTime = 0.0

  !! hand myPE out to C routines to avoid architecture-dependent code
  call driver_abortflashc_set_mype(dr_globalMe)

  !! make sure our stack (and whatever other rlimits) are big enough.  
  !! this should get around the 2Mb stack limit that pthreads
  !! imposes if linked statically (but not dynamically!)
  call dr_set_rlimits(dr_globalMe)

    
  !! Initialize runtime parameters
  call RuntimeParameters_init(dr_restart)

  call Driver_setupParallelEnv()

  !! Initialize the code timers.  Ideally should be first thing in
  !! code but currently the timing package
  !! uses MPI_WTime(), so Driver_initParallel() must go first, and
  !! uses RuntimeParameters_get(), so RuntimeParameters_init() must go
  !! first.
  call Profiler_init()
  call Timers_init(dr_initialWCTime)
  call Timers_start("initialization")

  !! Initialize the numerical tools.
  !! Right now this can go anywhere (NF)
  call Driver_initNumericalTools ()

  !PhysicalConstants init and Multispecies init must come before Logfile
  !since their values are stamped to the logfile
  call PhysicalConstants_init()

  !must come before EOS
  call Multispecies_init()

  call Logfile_init()

!!$  call Grid_init()
  
  call Driver_initMaterialProperties()
  if(dr_globalMe==MASTER_PE)print*,'MaterialProperties initialized'
  

  call RuntimeParameters_get('dtInit',dr_dtInit)

     call Eos_init()


  !Hydro_init must go before Driver
  if(dr_globalMe==MASTER_PE) print *, 'Ready to call Hydro_init' 
  call Hydro_init()           ! Hydrodynamics, MHD, RHD
  if(dr_globalMe==MASTER_PE)print*,'Hydro initialized'
  
 
  call Particles_initData(dr_restart,dr_particlesInitialized)
  

  !!Done with initialization.
  call Timers_stop ("initialization")

  call Eos_logDiagnostics(.TRUE.)

  return
end subroutine Driver_initFlash
