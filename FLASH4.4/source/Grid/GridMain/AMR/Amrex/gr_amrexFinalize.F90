!!****if* source/Grid/GridMain/AMR/Amrex/gr_amrexFinalize
!!
!! NAME
!!  gr_amrexFinalize
!!
!! SYNOPSIS
!!  gr_amrexFinalize
!!           
!! DESCRIPTION
!!  Clean-up all data structures managed by the Grid unit and allow AMReX to
!!  clean-up and terminate octree-based operations.
!!
!!***

subroutine gr_amrexFinalize()
    use iso_c_binding
    use amrex_init_module,      ONLY : amrex_finalize
    use amrex_amrcore_module,   ONLY : amrex_max_level, &
                                       amrex_amrcore_finalize
    use amrex_multifab_module,  ONLY : amrex_multifab_destroy
    use amrex_octree_module,    ONLY : amrex_octree_finalize

    use gr_physicalMultifabs,   ONLY : unk, &
                                       facevarx, facevary, facevarz
 
    integer :: lev

    write(*,*) "[gr_amrexFinalize] Finalizing"
  
    ! NOTE: Arrays of multifabs use AMReX's 0-based level indexing scheme
    do lev = 0, amrex_max_level
        call amrex_multifab_destroy(unk(lev))
        call amrex_multifab_destroy(facevarx(lev))
        call amrex_multifab_destroy(facevary(lev))
        call amrex_multifab_destroy(facevarz(lev))
    end do

    if (allocated(unk))         deallocate(unk)
    if (allocated(facevarx))    deallocate(facevarx)
    if (allocated(facevary))    deallocate(facevary)
    if (allocated(facevarz))    deallocate(facevarz)

    call amrex_amrcore_finalize()
    call amrex_octree_finalize()
    call amrex_finalize()
end subroutine gr_amrexFinalize

