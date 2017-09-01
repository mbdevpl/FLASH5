subroutine gr_clearLevel(lev) bind(c)
    use amrex_amr_module, ONLY : amrex_multifab_destroy
    use physicaldata,     ONLY : unk

    implicit none

    integer, intent(in), value :: lev

    ! C++ AMReX uses zero-based level index set, but FLASH uses 1-based set
    call amrex_multifab_destroy(unk(lev + 1))
end subroutine gr_clearLevel 

