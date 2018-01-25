!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsFinalize
!!
!! NAME
!!
!!  gr_amrexLsFinalize
!!
!!
!! SYNOPSIS
!!
!!  call gr_amrexLsFinalize ()
!!
!! Description
!!  
!!  Release HYPRE objects and other temporary storage from memory.
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!***

subroutine gr_amrexLsFinalize()
    use gr_amrexLsData, ONLY : gr_amrexLs_max_level, gr_amrexLs_geom,           &
                           gr_amrexLs_ba, gr_amrexLs_dm, gr_amrexLs_rhs, & 
                           gr_amrexLs_solution, gr_amrexLs_exact_solution
    use amrex_geometry_module, ONLY : amrex_geometry_destroy
    use amrex_boxarray_module, ONLY : amrex_boxarray_destroy
    use amrex_distromap_module, ONLY : amrex_distromap_destroy
    use amrex_multifab_module, ONLY : amrex_multifab_destroy

  implicit none
  !!Nothing here yet to finalize
    integer :: ilev, idim
    do ilev = 0, gr_amrexLs_max_level
       call amrex_geometry_destroy(gr_amrexLs_geom(ilev))
       call amrex_boxarray_destroy(gr_amrexLs_ba(ilev))
       call amrex_distromap_destroy(gr_amrexLs_dm(ilev))
       call amrex_multifab_destroy(gr_amrexLs_solution(ilev))
       call amrex_multifab_destroy(gr_amrexLs_rhs(ilev))
       call amrex_multifab_destroy(gr_amrexLs_exact_solution(ilev))
!        !Abec solver not implemented yet
!        if (allocated(acoef)) then
!           call amrex_multifab_destroy(acoef(ilev))
!           call amrex_multifab_destroy(bcoef(ilev))
!        end if
    end do
  
end subroutine gr_amrexLsFinalize