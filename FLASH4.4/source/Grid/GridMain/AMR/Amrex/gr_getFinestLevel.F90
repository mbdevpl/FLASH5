subroutine gr_getFinestLevel(level)
    use amrex_amrcore_module, ONLY : amrex_get_finest_level

    implicit none

    integer, intent(OUT) :: level

    ! AMReX uses 0-based level indexing/FLASH uses 1-based
    level = amrex_get_finest_level() + 1
end subroutine gr_getFinestLevel

