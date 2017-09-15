subroutine gr_remakeLevel(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,      ONLY : wp => amrex_real
    use amrex_box_module,       ONLY : amrex_box
    use amrex_boxarray_module,  ONLY : amrex_boxarray, &
                                       box_print => amrex_print
    use amrex_distromap_module, ONLY : amrex_distromap, &
                                       distro_print => amrex_print
    use amrex_parallel_module,  ONLY : amrex_parallel_myproc
    use amrex_interfaces,       ONLY : gr_clearLevel

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm
 
    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_box)       :: bx

    integer :: rank = 0

    ba = pba
    dm = pdm
 
    rank = amrex_parallel_myproc()
    write(*,*) "[Rank ", rank, "] gr_remakeLevel - Start Level ", lev + 1

    call gr_clearLevel(lev)

    write(*,*) "[Rank ", rank, "] gr_remakeLevel - End Level ", lev + 1
end subroutine gr_remakeLevel 

