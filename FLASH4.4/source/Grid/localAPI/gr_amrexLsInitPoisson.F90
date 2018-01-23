!!****if* source/Grid/GridSolvers/HYPRE/gr_amrexLsInitPoisson
!!
!!  NAME 
!!
!! gr_amrexLsInitPoisson
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsInitPoisson()
!!
!!
!!  DESCRIPTION 
!! This routine sets up the Poisson solver from the 
!! Amrex Linear Solvers. This is called from  Grid_solvePoisson
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

subroutine gr_amrexLsInitPoisson (geom, solution, rhs, exact_solution)
  
        use amrex_geometry_module, ONLY : amrex_geometry
        use amrex_multifab_module, ONLY : amrex_multifab
        implicit none
        type(amrex_geometry), intent(in   ) :: geom(0:)
        type(amrex_multifab), intent(inout) :: solution(0:)
        type(amrex_multifab), intent(inout) :: rhs(0:)
        type(amrex_multifab), intent(inout) :: exact_solution(0:)

   
  return 
  
end subroutine gr_amrexLsInitPoisson
