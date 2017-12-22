subroutine gr_amrextAmrexFinalize()
    use iso_c_binding
    use amrex_amr_module
    use amrex_octree_module,    ONLY : amrex_octree_finalize
    use gr_physicalMultifabs,   ONLY : unk, &
                                       facevarx, facevary, facevarz
    use gr_amrextData,          ONLY : gr_amrextDataFree
 
    integer :: lev

    write(*,*) "[gr_amrextAmrexFinalize] Finalizing"
  
    ! NOTE: We implement these with the 1-based level indexing scheme native to
    ! so that the AMReX unk has a similar interface to the paramesh unk.
    !   => all code dealing with multifabs at the Fortran/C++ interface must take 
    !      care of the index translation
    do lev = 0, amrex_max_level
       if (lev .le. ubound(unk,1)) call amrex_multifab_destroy(unk(lev))
    end do

    deallocate(unk)
    if (allocated(facevarx)) deallocate(facevarx)
    if (allocated(facevary)) deallocate(facevary)
    if (allocated(facevarz)) deallocate(facevarz)

    call amrex_amrcore_finalize()
    call amrex_octree_finalize()
    call gr_amrextDataFree()
    if (amrex_initialized()) call amrex_finalize()
end subroutine gr_amrextAmrexFinalize

