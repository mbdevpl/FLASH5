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
    use amrex_fillpatch_module,    ONLY : amrex_fillpatch
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons

    use Grid_data,                 ONLY : lo_bc_amrex, hi_bc_amrex
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

    integer :: j

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
       call amrex_fillpatch(mfab, time+1.0d0, unk(lev), &
                                  time,       unk(lev), &
                                  amrex_geom(lev), gr_fillPhysicalBC, &
                                  time, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS)       
    else
       call amrex_fillpatch(mfab, time+1.0d0, unk(lev-1), &
                                  time,       unk(lev-1), &
                                  amrex_geom(lev-1), gr_fillPhysicalBC, &
                                  time+1.0e0, unk(lev  ), &
                                  time,       unk(lev  ), &
                                  amrex_geom(lev  ), gr_fillPhysicalBC, &
                                  time, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                  amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                  lo_bc_amrex, hi_bc_amrex)       
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

