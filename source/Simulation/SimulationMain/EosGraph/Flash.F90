!!****if* source/Simulation/SimulationMain/EosGraph/Flash
!!
!! NAME
!!
!!  Flash
!!
!!
!! SYNOPSIS
!!
!!  N/A
!!
!!
!! DESCRIPTION
!!
!!  The source file Flash.F90 in the Simulation unit contains the Fortran
!!  PROGRAM. As such it can be considered the top-level "driver" of an application.
!!  By default it is set up to drive the simulation of a time-dependent
!!  problem by calling:
!!  - Driver_initFlash  for initializations,
!!  - Driver_evolveFlash  for managing the computation, and
!!  - Driver_finalizeFlash  for cleaning up.
!!
!! SEE ALSO
!!
!!  Driver_initFlash
!!  Driver_evolveFlash
!!  Driver_finalizeFlash
!!
!!***

program Flash

  use Driver_interface, ONLY : Driver_initParallel, Driver_initFlash,&
       Driver_evolveFlash

  implicit none

  call Driver_initParallel()

  call Driver_initFlash()

  call Driver_evolveFlash( )


end program Flash
