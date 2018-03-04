#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_clearLevelCallback(lev) bind(c)
    use amrex_amr_module,          ONLY : amrex_multifab_destroy
    use amrex_fluxregister_module, ONLY : amrex_fluxregister_destroy 

    use Grid_data,                 ONLY : gr_doFluxCorrection
    use gr_physicalMultifabs,      ONLY : unk, &
                                          facevarx, facevary, facevarz, &
                                          fluxes, &
                                          flux_registers

    implicit none

    integer, intent(in), value :: lev

    integer :: dir

    ! Multifab arrays use 0-based index set like AMReX
    call amrex_multifab_destroy(unk     (lev))
#if NFACE_VARS > 0
    call amrex_multifab_destroy(facevarx(lev))
    call amrex_multifab_destroy(facevary(lev))
    call amrex_multifab_destroy(facevarz(lev))
#endif

#if NFLUXES > 0
    do dir = 1, SIZE(fluxes, 2)
        call amrex_multifab_destroy(fluxes(lev, dir))
    end do

    if ((lev > 0) .AND. (gr_doFluxCorrection)) then
        call amrex_fluxregister_destroy(flux_registers(lev))
    end if
#endif

#ifdef DEBUG_GRID
    write(*,'(A,A,I2)') "[gr_clearLevelCallback]", &
                        "              Cleared level", lev + 1
#endif

end subroutine gr_clearLevelCallback

