subroutine gr_amrex_finalize()
    use iso_c_binding
    use amrex_amr_module
    use amrex_octree_module,    ONLY : amrex_octree_finalize
    use physicaldata,           ONLY : unk
 
    integer :: lev

    write(*,*) "[gr_amrex_finalize] Finalizing"
  
    ! NOTE: We implement these with the 1-based level indexing scheme native to
    ! so that the AMReX unk has a similar interface to the paramesh unk.
    !   => all code dealing with multifabs at the Fortran/C++ interface must take 
    !      care of the index translation
    do lev = 1, (amrex_max_level + 1)
        call amrex_multifab_destroy(unk(lev))
    end do

    deallocate(unk)

    call amrex_amrcore_finalize()
    call amrex_octree_finalize()
    call amrex_finalize()
end subroutine gr_amrex_finalize

