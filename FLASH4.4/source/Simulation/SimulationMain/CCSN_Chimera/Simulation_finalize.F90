!!****f* source/Simulation/Simulation_finalize
!!
!! NAME
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This function cleans up the Simulation unit, deallocates memory, etc.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Simulation_finalize()

  use chimera_model_module

  implicit none

  call close_chimera_file

end subroutine Simulation_finalize
