subroutine gr_clearLevel(lev) bind(c)
    use amrex_amr_module, ONLY : amrex_multifab_destroy
    use physicaldata,     ONLY : unk

    implicit none

    integer, intent(in), value :: lev

    call amrex_multifab_destroy(unk(lev))
end subroutine gr_clearLevel 

