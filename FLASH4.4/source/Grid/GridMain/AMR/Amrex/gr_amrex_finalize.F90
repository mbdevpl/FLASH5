subroutine gr_amrex_finalize()
    use iso_c_binding
    use amrex_amr_module
    use amrex_octree_module,    ONLY : amrex_octree_finalize
    use physicaldata,           ONLY : unk
 
    integer :: lev
    do lev = 0, amrex_max_level
        call amrex_multifab_destroy(unk(lev))
    end do

    deallocate(unk)

    call amrex_amrcore_finalize()
    call amrex_octree_finalize()
    call amrex_finalize()
end subroutine gr_amrex_finalize

