#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Flash.h"
#include "constants.h"

subroutine gr_remakeLevelCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio, &
                                          amrex_geom
    use amrex_box_module,          ONLY : amrex_box
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_multifab, &
                                          amrex_multifab_build, &
                                          amrex_multifab_destroy
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons

    use Grid_data,                 ONLY : gr_lo_bc_ptr, gr_hi_bc_ptr
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
    type(amrex_box)       :: bx
    type(amrex_multifab)  :: mfab

    integer, target :: lo_bc(NDIM, UNK_VARS_BEGIN:UNK_VARS_END)
    integer, target :: hi_bc(NDIM, UNK_VARS_BEGIN:UNK_VARS_END)
    type(c_ptr)     :: lo_bc_ptr(UNK_VARS_BEGIN:UNK_VARS_END)
    type(c_ptr)     :: hi_bc_ptr(UNK_VARS_BEGIN:UNK_VARS_END)

    type(c_ptr) :: mfab_src(1)
    real(wp)    :: time_src(1)
    type(c_ptr) :: mfab_coarse(1)
    real(wp)    :: time_coarse(1)
    type(c_ptr) :: mfab_fine(1)
    real(wp)    :: time_fine(1)

    integer :: j

    ! AMReX C++ fillpatch routines
    interface
       subroutine amrex_fi_fillpatch_single(mf, time, smf, stime, ns, scomp, dcomp, ncomp, &
                                            geom, fill) bind(c)
         import
         implicit none
         type(c_ptr),      value :: mf
         real(wp),         value :: time
         type(c_ptr), intent(in) :: smf(*)
         real(wp),    intent(in) :: stime(*)
         integer(c_int),   value :: scomp, dcomp, ncomp, ns
         type(c_ptr),      value :: geom
         type(c_funptr),   value :: fill
       end subroutine amrex_fi_fillpatch_single

       subroutine amrex_fi_fillpatch_two(mf, time, &
                                         cmf, ctime, nc, fmf, ftime, nf, scomp, dcomp, ncomp, &
                                         cgeom, fgeom, cfill, ffill, rr, interp, lo_bc, hi_bc) bind(c)
         import
         implicit none
         type(c_ptr),      value :: mf
         real(wp),         value :: time
         type(c_ptr), intent(in) :: cmf(*)
         real(wp),    intent(in) :: ctime(*)
         integer,          value :: nc
         type(c_ptr), intent(in) :: fmf(*)
         real(wp),    intent(in) :: ftime(*)
         integer,          value :: nf
         integer,          value :: scomp, dcomp, ncomp
         type(c_ptr),      value :: cgeom, fgeom
         type(c_funptr),   value :: cfill, ffill
         integer,          value :: rr, interp
         type(c_ptr), intent(in) :: lo_bc(*), hi_bc(*)
       end subroutine amrex_fi_fillpatch_two
    end interface

    ba = pba
    dm = pdm

#ifdef DEBUG_GRID
    write(*,'(A,A,I2)') "[gr_remakeLevelCallback]", &
                      "             Start Level ", lev + 1
#endif

    !!!!! SAVE DATA IN BUFFER WITH GIVEN BOXARRAY/DISTRIBUTION
    ! Get all unk interior data
    call amrex_multifab_build(mfab, ba, dm, NUNK_VARS, NGUARD)
    ! DEVNOTE: TODO Include facevars in this process

    if (lev == 0) then
       ! Move all unk data to given ba/dm layout.  Do *not* use sub-cycling.
       ! -1 because of Fortran variable index starts with 1
       mfab_src(1) = unk(lev)%p
       time_src(1) = time
       call amrex_fi_fillpatch_single(mfab%p, time, mfab_src, time_src, 1, &
                                      UNK_VARS_BEGIN-1, UNK_VARS_BEGIN-1, &
                                      NUNK_VARS, amrex_geom(lev)%p, &
                                      c_funloc(gr_fillPhysicalBC))
    else
       mfab_coarse(1) = unk(lev-1)%p
       time_coarse(1) = time
       mfab_fine(1)   = unk(lev  )%p
       time_fine(1)   = time

       call amrex_fi_fillpatch_two(mfab%p, time, &
                                   mfab_coarse, time_coarse, 1, &
                                   mfab_fine, time_fine, 1, &
                                   UNK_VARS_BEGIN-1, UNK_VARS_BEGIN-1, NUNK_VARS, &
                                   amrex_geom(lev-1)%p, amrex_geom(lev)%p, &
                                   c_funloc(gr_fillPhysicalBC), &
                                   c_funloc(gr_fillPhysicalBC), &
                                   amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                   gr_lo_bc_ptr, gr_hi_bc_ptr)
    end if

    !!!!! REBUILD MFAB AT LEVEL AND FILL FROM BUFFER
    call gr_clearLevelCallback(lev)
    call amrex_multifab_build(unk(lev), ba, dm, NUNK_VARS, NGUARD)

    ! Only copy interior
    call unk(lev)%copy(mfab, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                       NGUARD)

    call amrex_multifab_destroy(mfab)

    write(*,'(A,A,I2)') "[gr_remakeLevelCallback]", &
                      "             Remake level ", lev + 1

end subroutine gr_remakeLevelCallback

