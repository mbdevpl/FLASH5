!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsSolvePoissonUnk
!!
!!  NAME 
!!
!! gr_amrexLsSolvePoissonUnk
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsSolvePoissonUnk()
!!
!!
!!  DESCRIPTION 
!! This routine solves the Poisson equation from the 
!! Amrex Linear Solvers using the variables from Unk multifab
!! for rhs and unknown phi 
!!
!!
!! ARGUMENTS
!! 
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Solver settings used ?? DESCRIPTION??
!!
!!***

subroutine gr_amrexLsSolvePoissonUnk ()
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use amrex_multigrid_module, ONLY : amrex_multigrid, amrex_multigrid_build, amrex_multigrid_destroy
  use amrex_poisson_module, ONLY : amrex_poisson, amrex_poisson_build, amrex_poisson_destroy
  use amrex_lo_bctypes_module, ONLY : amrex_lo_periodic

  implicit none
  
#include "Flash.h"
#include "constants.h"   
  
  call Timers_start("gr_amrexLsInitPoissonUnk")
  

  call Timers_stop("gr_amrexLsInitPoissonUnk")
end subroutine gr_amrexLsSolvePoissonUnk
