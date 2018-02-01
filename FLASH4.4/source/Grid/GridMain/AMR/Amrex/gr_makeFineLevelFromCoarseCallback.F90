#include "constants.h"
#include "Flash.h"

subroutine gr_makeFineLevelFromCoarseCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio, &
                                          amrex_geom
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_multifab_build, &
                                          amrex_mfiter, &
                                          amrex_mfiter_build, &
                                          amrex_mfiter_destroy
    use amrex_fillpatch_module,    ONLY : amrex_fillcoarsepatch
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons

    use Grid_data,                 ONLY : lo_bc_amrex, hi_bc_amrex
    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_fillPhysicalBC
    use gr_physicalMultifabs,      ONLY : unk, &
                                          facevarx, facevary, facevarz
    use Driver_interface,          ONLY : Driver_abortFlash

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter)    :: mfi

    integer :: nFab

    call Driver_abortFlash("[gr_makeFileLevelFromCoarseCallback] " // &
                           "Callback has never been tested")

    ba = pba
    dm = pdm

    !!!!!----- (Re)create FABS for storing physical data at this level
    call gr_clearLevelCallback(lev)
    call amrex_multifab_build(unk     (lev), ba, dm, NUNK_VARS, NGUARD)
    ! DEVNOTE: TODO Create these wrt proper face-centered boxes
    call amrex_multifab_build(facevarx(lev), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevary(lev), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevarz(lev), ba, dm, NUNK_VARS, NGUARD)

    !!!!!----- Fill new refinement level via interpolation from parent block
    ! This *hopefully* will do the guard cell fill as well
    ! NOTE: FLASH does not use sub-cycling (temporal interpolation)
    !
    ! -1 because Fortran variable index starts with 1
    call amrex_fillcoarsepatch(unk(lev), time,     unk(lev-1),  &
                                         time+0.1, unk(lev-1),  &
                                         amrex_geom(lev-1), gr_fillPhysicalBC,  &
                                         amrex_geom(lev  ), gr_fillPhysicalBC,  &
                                         time, &
                                         UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                         amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                         lo_bc_amrex, hi_bc_amrex) 

    ! DEV: FIXME Should we do an EoS run on interiors and GC here?

    nFab = 0
    call amrex_mfiter_build(mfi, unk(lev), tiling=.false.)
    do while(mfi%next())
        nFab = nFab + 1 
    end do
    call amrex_mfiter_destroy(mfi)

    write(*,'(A,I0,A,I0,A)') "Made fine level ", lev + 1, " - ", nFab, " blocks"

end subroutine gr_makeFineLevelFromCoarseCallback

