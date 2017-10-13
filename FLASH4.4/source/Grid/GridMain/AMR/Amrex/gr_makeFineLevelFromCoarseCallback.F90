#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine gr_makeFineLevelFromCoarseCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio, &
                                          amrex_geom
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_multifab_build
    use amrex_bc_types_module,     ONLY : amrex_bc_int_dir
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons

    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_fillPhysicalBC
    use gr_physicalMultifabs,      ONLY : unk, &
                                          facevarx, facevary, facevarz

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm

    integer, target :: lo_bc(NDIM, UNK_VARS_BEGIN:UNK_VARS_END)
    integer, target :: hi_bc(NDIM, UNK_VARS_BEGIN:UNK_VARS_END)
    type(c_ptr)     :: lo_bc_ptr(UNK_VARS_BEGIN:UNK_VARS_END)
    type(c_ptr)     :: hi_bc_ptr(UNK_VARS_BEGIN:UNK_VARS_END)

    integer :: j

    ! AMReX C++ fillpatch routines
    interface
      subroutine amrex_fi_fillcoarsepatch(mf, time, cmf, scomp, dcomp, ncomp, &
                                         cgeom, fgeom, cfill, ffill, rr, &
                                         interp, lo_bc, hi_bc) bind(c)
         import
         implicit none
         type(c_ptr), value :: mf, cmf, cgeom, fgeom
         type(c_ptr), intent(in) :: lo_bc(*), hi_bc(*)
         type(c_funptr), value :: cfill, ffill
         real(wp), value :: time
         integer, value :: scomp, dcomp, ncomp, rr, interp
      end subroutine amrex_fi_fillcoarsepatch
    end interface

#ifdef DEBUG_GRID
    write(*,'(A,I2)') "[gr_makeFineLevelFromCoarseCallback] Start on level ", lev + 1
#endif

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
    ! DEVNOTE: FIXME Currently fixing BC to periodic here
    ! DEVNOTE: FIXME Currently fixing interpolation mode to cell conserved
    !                linear (AMReX_Interpolater.H)
    lo_bc(:, :) = amrex_bc_int_dir
    hi_bc(:, :) = amrex_bc_int_dir
    do j = UNK_VARS_BEGIN, UNK_VARS_END
       lo_bc_ptr(j) = c_loc(lo_bc(1, j))
       hi_bc_ptr(j) = c_loc(hi_bc(1, j))
    end do

    ! -1 because Fortran variable index starts with 1
    call amrex_fi_fillcoarsepatch(unk(lev)%p, time, unk(lev-1)%p, &
                                  UNK_VARS_BEGIN-1, UNK_VARS_BEGIN-1, NUNK_VARS, &
                                  amrex_geom(lev-1)%p, amrex_geom(lev)%p, &
                                  c_funloc(gr_fillPhysicalBC), &
                                  c_funloc(gr_fillPhysicalBC), &
                                  amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                  lo_bc_ptr, hi_bc_ptr)

    write(*,'(A,I2)') "[gr_makeFineLevelFromCoarseCallback] Make fine level ", lev + 1

end subroutine gr_makeFineLevelFromCoarseCallback

