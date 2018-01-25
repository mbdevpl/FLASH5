!!****ih* source/Grid/GridSolvers/AmrexMultiGrid/gr_amrexLsInterface
!!
!! This is a module file for the AmrexMultiGrid solver in FLASH that defines
!! additional interfaces private to the GridSolvers/AmrexMultiGrid implementations.
!!
!! NOTES
!!
!!  Adding explicit interfaces here is being tried as an alternative to
!!  writiting executable FORTRAN wrappers for each additional HYPRE routine
!!  we want to call from FLASH.
!!
!!***


Module gr_amrexLsInterface
!     use amrex_fort_module,     ONLY : amrex_real
!     use amrex_box_module,     ONLY : amrex_box
!     use amrex_boxarray_module,     ONLY : amrex_boxarray
!     use amrex_distromap_module,     ONLY : amrex_distromap
!     use amrex_geometry_module, ONLY : amrex_geometry
!     use amrex_multifab_module, ONLY : amrex_multifab, amrex_mfiter

    implicit none
    
#include "constants.h"
    !
    interface
    subroutine gr_amrexLsInitPoisson (geom, solution, rhs, exact_solution)
        use amrex_geometry_module, ONLY : amrex_geometry
        use amrex_multifab_module, ONLY : amrex_multifab
        implicit none
        type(amrex_geometry), intent(in   ) :: geom(0:)
        type(amrex_multifab), intent(inout) :: solution(0:)
        type(amrex_multifab), intent(inout) :: rhs(0:)
        type(amrex_multifab), intent(inout) :: exact_solution(0:)
    end subroutine gr_amrexLsInitPoisson
    end interface

    interface
    subroutine gr_amrexLsSolvePoisson ()
        implicit none
    end subroutine gr_amrexLsSolvePoisson
    end interface

    interface
    subroutine gr_amrexLsInitGeom ()
        implicit none
    end subroutine gr_amrexLsInitGeom
    end interface

  interface
    subroutine gr_amrexLsInitGrid ()
        implicit none
    end subroutine gr_amrexLsInitGrid
  end interface

  interface
    subroutine gr_amrexLsInitMf ()
        implicit none
    end subroutine gr_amrexLsInitMf
  end interface
!   interface
!      subroutine gr_hypreGetFaceBFcB (direction, blkLimits, blkLimitsGC, facBptr, flux, iVar)
!        implicit none
!        integer, intent(IN) :: direction
!        integer, intent(IN) :: blkLimits (2,MDIM) 
!        integer, intent(IN) :: blkLimitsGC (2,MDIM)
!        real, POINTER_INTENT_IN :: facBptr(:,:,:)   
!        real, intent(INOUT) :: flux(:,:,:,:)
!        integer, intent(IN) :: iVar
!      end subroutine gr_hypreGetFaceBFcB
!   end interface

end Module gr_amrexLsInterface
