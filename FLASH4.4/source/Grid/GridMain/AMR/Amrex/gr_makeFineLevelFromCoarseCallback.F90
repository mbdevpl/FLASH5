subroutine gr_makeFineLevelFromCoarseCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio, &
                                          amrex_geom
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_fillpatch_module,    ONLY : amrex_fillcoarsepatch
    use amrex_multifab_module,     ONLY : amrex_multifab_build
    use amrex_bc_types_module,     ONLY : amrex_bc_int_dir
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons

    use amrex_interfaces,          ONLY : gr_clearLevelCallback, &
                                          gr_fillPhysicalBC
    use Grid_data,                 ONLY : gr_iguard
    use gr_physicalMultifabs,      ONLY : unk, &
                                          facevarx, facevary, facevarz

    implicit none

#include "constants.h"
#include "Flash.h"

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm

    integer :: lo_bc(NDIM, 1)
    integer :: hi_bc(NDIM, 1)

    write(*,*) "[gr_makeFineLevelFromCoarseCallback] Start on level ", lev + 1

    ba = pba
    dm = pdm

    !!!!!----- (Re)create FABS for storing physical data at this level
    call gr_clearLevelCallback(lev)
    call amrex_multifab_build(unk     (lev), ba, dm, NUNK_VARS, gr_iguard)
    ! DEVNOTE: TODO Create these wrt proper face-centered boxes
    call amrex_multifab_build(facevarx(lev), ba, dm, NUNK_VARS, gr_iguard)
    call amrex_multifab_build(facevary(lev), ba, dm, NUNK_VARS, gr_iguard)
    call amrex_multifab_build(facevarz(lev), ba, dm, NUNK_VARS, gr_iguard)

    !!!!!----- Fill new refinement level via interpolation from parent block
    ! This *hopefully* will do the guard cell fill as well
    ! NOTE: FLASH does not use sub-cycling (temporal interpolation)
    !
    ! DEVNOTE: FIXME I think the BC here should be arrays
    ! DEVNOTE: FIXME Currently fixing BC to periodic here
    ! DEVNOTE: FIXME Currently fixing interpolation mode to cell conserved
    !                linear (AMReX_Interpolater.H)
    ! DEVNOTE: TODO Since we are not using subcycling, should we just use
    !               amrex_fi_fillcoarsepatch directly?
    lo_bc(:, :) = amrex_bc_int_dir
    hi_bc(:, :) = amrex_bc_int_dir
    call amrex_fillcoarsepatch(unk(lev), time,     unk(lev-1),  &
                                         time+0.1, unk(lev-1),  &
                               amrex_geom(lev-1), gr_fillPhysicalBC,  &
                               amrex_geom(lev  ), gr_fillPhysicalBC,  &
                               time, &
                               UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                               amrex_ref_ratio(lev), amrex_interp_cell_cons, &
                               lo_bc, hi_bc)

    write(*,*) "[gr_makeFineLevelFromCoarseCallback] Finished on level ", lev + 1
end subroutine gr_makeFineLevelFromCoarseCallback

