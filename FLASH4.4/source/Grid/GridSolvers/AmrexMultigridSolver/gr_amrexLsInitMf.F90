!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsInitMf
!!
!!  NAME 
!!
!! gr_amrexLsInitMf
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsInitMf()
!!
!!
!!  DESCRIPTION 
!! This routine sets up the Multifab from the 
!! Amrex Linear Solvers. 
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

subroutine gr_amrexLsInitMf ()
  
!   use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_amrexLsData, ONLY : gr_amrexLs_ba, gr_amrexLs_dm, &
!                                             gr_amrexLs_ascalar, gr_amrexLs_bscalar, &
                                            gr_amrexLs_max_level, gr_amrexLs_prob_type, &
                                            gr_amrexLs_exact_solution, gr_amrexLs_rhs, gr_amrexLs_solution
  use amrex_fort_module,     ONLY : amrex_real
  use amrex_distromap_module,     ONLY : amrex_distromap_build
  use amrex_multifab_module,    ONLY : amrex_multifab_build
  implicit none
  
#include "Flash.h"
#include "constants.h"   
  
    integer :: ilev

    call Timers_start("gr_amrexLsInitMf")     
  
    do ilev = 0, gr_amrexLs_max_level
       call amrex_distromap_build(gr_amrexLs_dm(ilev),gr_amrexLs_ba(ilev))
       ! one ghost cell to store boundary conditions
       call amrex_multifab_build(gr_amrexLs_solution(ilev), gr_amrexLs_ba(ilev), gr_amrexLs_dm(ilev), nc=1, ng=1)
       call amrex_multifab_build(gr_amrexLs_rhs(ilev), gr_amrexLs_ba(ilev), gr_amrexLs_dm(ilev), nc=1, ng=0)
       call amrex_multifab_build(gr_amrexLs_exact_solution(ilev), gr_amrexLs_ba(ilev), gr_amrexLs_dm(ilev), nc=1, ng=0)
!      For abec. Not yet implemented
!        if (gr_amrexLs_prob_type==2 && (gr_amrexLs_acoef)) then
!           call amrex_multifab_build(gr_amrexLs_acoef(ilev), gr_amrexLs_ba(ilev), gr_amrexLs_dm(ilev), nc=1, ng=0)
!           ! 1 ghost cell for averaging from cell centers to faces
!           call amrex_multifab_build(gr_amrexLs_bcoef(ilev), gr_amrexLs_ba(ilev), gr_amrexLs_dm(ilev), nc=1, ng=1)
!        end if
    end do
  
  call Timers_stop("gr_amrexLsInitMf")    

end subroutine gr_amrexLsInitMf
