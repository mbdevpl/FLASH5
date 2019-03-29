!!****if* source/Driver/DriverMain/Driver_finalizeFlash
!!
!! NAME
!!  Driver_finalizeFlash
!!
!! SYNOPSIS
!!  Driver_finalizeFlash()
!!
!! DESCRIPTION
!!
!!  Calls all the unit finalize routines
!!  which may need
!!  memory deallocated etc before the run end.
!!  Order does matter.
!!
!!***


subroutine Driver_finalizeFlash()

  use Eos_interface, ONLY : Eos_finalize
  use Burn_interface, ONLY : Burn_finalize
  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Multispecies_interface, ONLY : Multispecies_finalize
  use Particles_interface, ONLY : Particles_finalize
  use Grid_interface, ONLY : Grid_finalize
  use Hydro_interface, ONLY : Hydro_finalize
  use Driver_data, ONLY: dr_globalMe, dr_restart
  use Simulation_interface, ONLY : Simulation_finalize
  use IO_interface, ONLY : IO_finalize
  use Timers_interface, ONLY: Timers_finalize
  use Gravity_interface, ONLY: Gravity_finalize
implicit none
#include "mpif.h"

  integer :: ierr
 
  
!!$  call Profiler_finalize()
!!$  
  call RuntimeParameters_finalize()

  call Multispecies_finalize()

  call Burn_finalize()             ! NSE material property

  call Grid_finalize()            ! Grid package
 
  call Particles_finalize()       ! Particles
  
  call Hydro_finalize()           ! Hydrodynamics
  
  call Eos_finalize()             ! Equation of State

  call Gravity_finalize()         ! Gravity
  call IO_finalize()

  call Simulation_finalize()

  call Timers_finalize()

  call MPI_Finalize(ierr)

  return
end subroutine Driver_finalizeFlash








