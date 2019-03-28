!!****if* source/Grid/GridMain/AMR/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  A stub implementation as this public Grid routine is not needed for
!!  AMReX-based FLASH simulations.
!!
!!  Refinement of blocks is managed by AMReX through through the
!!  gr_markRefineDerefineCallback routine.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Grid_markRefineDerefine()
  use Driver_interface, ONLY : Driver_abortFlash
  
  implicit none
  
  call Driver_abortFlash("[Grid_markRefineDerefine] Not implemented in AMReX")
end subroutine Grid_markRefineDerefine

