subroutine gr_clearLevel(lev) bind(c)
    use amrex_amr_module, ONLY : amrex_multifab_destroy
    use gr_physicalMultifabs,   ONLY : unk, &
                                       facevarx, facevary, facevarz

    implicit none

    integer, intent(in), value :: lev
  
    ! Multifab arrays use 0-based index set like AMReX
    write(*,*) "[gr_clearLevel] called on level", lev + 1
    call amrex_multifab_destroy(unk     (lev))
    call amrex_multifab_destroy(facevarx(lev))
    call amrex_multifab_destroy(facevary(lev))
    call amrex_multifab_destroy(facevarz(lev))
end subroutine gr_clearLevel 

