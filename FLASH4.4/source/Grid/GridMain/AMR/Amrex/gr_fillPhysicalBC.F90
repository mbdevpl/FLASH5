#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_fillPhysicalBC(pmf, scomp, ncomp, time, pgeom) bind(c)
    use iso_c_binding
    use amrex_fort_module, ONLY : wp => amrex_real

    implicit none

    type(c_ptr),    value :: pmf
    type(c_ptr),    value :: pgeom
    integer(c_int), value :: scomp
    integer(c_int), value :: ncomp
    real(wp),       value :: time

    ! In this test problem, we only have periodic boundaries.
    ! So there is nothing to do.
#ifdef DEBUG_GRID
    write(*,'(A,A,F7.4)') "[gr_fillPhysicalBC]", &
                        "                  Start for time ", time
    write(*,'(A,A,F7.4)') "[gr_fillPhysicalBC]", &
                        "                  Finished for time ", time
#endif

end subroutine gr_fillPhysicalBC

