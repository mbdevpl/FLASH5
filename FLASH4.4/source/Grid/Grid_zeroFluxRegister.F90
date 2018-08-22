!!****if* source/Grid/Grid_zeroFluxRegister
!!
!! NAME
!!  Grid_zeroFluxRegister
!!
!! SYNOPSIS
!!  call Grid_zeroFluxRegister(integer(IN) :: fine_level)
!!
!! DESCRIPTION 
!!  Each flux register is associated with a fine and a coarse level.  Given an
!!  index for the fine level, set the fine and coarse data to zero in the
!!  associated flux register.
!!
!! ARGUMENTS
!!  fine_level - the 1-based level index (1 is the coarsest level) that
!!               indicates the fine level of the flux register on which to 
!!               operate.
!!
!! SEE ALSO
!!   Grid_addFineToFluxRegister
!!   Grid_addCoarseToFluxRegister
!!
!!***

subroutine Grid_zeroFluxRegister(fine_level)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(IN) :: fine_level

  call Driver_abortFlash("[Grid_zeroFluxRegister] Prototype stub.  Do NOT use!")
end subroutine Grid_zeroFluxRegister

