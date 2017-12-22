#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine gr_initNewLevelCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amr_module,          ONLY : amrex_geom, &
                                          amrex_problo
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio
    use amrex_box_module,          ONLY : amrex_box
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_mfiter, &
                                          amrex_mfiter_build, &
                                          amrex_mfiter_destroy, &
                                          amrex_multifab_build
    use amrex_fillpatch_module,    ONLY : amrex_fillpatch
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons
    
    use gr_physicalMultifabs,      ONLY : unk, &
                                          facevarx, facevary, facevarz
    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_fillPhysicalBC
    use block_metadata,            ONLY : block_metadata_t
    use Simulation_interface,      ONLY : Simulation_initBlock
    use Grid_data,                 ONLY : gr_eosModeInit, &
                                          lo_bc_amrex, hi_bc_amrex
    use Eos_interface,             ONLY : Eos_wrapped

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter)    :: mfi
    type(amrex_box)       :: bx

    type(block_metadata_t)        :: block
    real(wp), contiguous, pointer :: initData(:,:,:,:)

    integer :: n_blocks

    ba = pba
    dm = pdm

    call gr_clearLevelCallback(lev)

    ! Create FABS for storing physical data at given level
    call amrex_multifab_build(unk     (lev), ba, dm, NUNK_VARS, NGUARD)
    ! DEVNOTE: TODO Create there w.r.t. proper face-centered boxes
#if NFACE_VARS > 0
    call amrex_multifab_build(facevarx(lev), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevary(lev), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevarz(lev), ba, dm, NUNK_VARS, NGUARD)
#endif

    ! Write initial data across domain at coarsest level
    call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)

    n_blocks = 0
    do while (mfi%next())
        bx = mfi%fabbox()

        ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
        ! Level must be 1-based index and limits/limitsGC must be 1-based also
        block%level = lev + 1
        block%grid_index = mfi%grid_index()
        block%limits(LOW,  :) = 1
        block%limits(HIGH, :) = 1
        block%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1 + NGUARD
        block%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1 - NGUARD
        block%limitsGC(LOW,  :) = 1
        block%limitsGC(HIGH, :) = 1
        block%limitsGC(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
        block%limitsGC(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1

        associate(lo => block%limitsGC(LOW, :))
            initData(lo(1):, lo(2):, lo(3):, 1:) => unk(lev)%dataptr(mfi)
        end associate
 
        !  We need to zero data in case we reuse blocks from previous levels
        !  but don't initialize all data in Simulation_initBlock... in particular
        !  the total vs. internal energies can cause problems in the eos call that 
        !  follows.
        initData = 0.0d0
        call Simulation_initBlock(initData, block)
        call Eos_wrapped(gr_eosModeInit, block%limits, initData)
        nullify(initData)

        n_blocks = n_blocks + 1
    end do

    call amrex_mfiter_destroy(mfi)

    ! Subsequent AMReX calls to gr_markRefineDerefineCallback require that the
    ! GC be filled.  We do *not* ask client code to do this, so fill GC here
    if (lev == 0) then
       ! Move all unk data to given ba/dm layout.  Do *not* use sub-cycling.
       ! -1 because of Fortran variable index starts with 1
       call amrex_fillpatch(unk(lev), time+1.0d0, unk(lev), &
                                      time,       unk(lev), &
                                      amrex_geom(lev), gr_fillPhysicalBC, &
                                      time, &
                                      UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS)
    else
       call amrex_fillpatch(unk(lev), time+1.0d0, unk(lev-1), &
                                      time,       unk(lev-1), &
                                      amrex_geom(lev-1), gr_fillPhysicalBC, &
                                      time+1.0e0, unk(lev  ), &
                                      time,       unk(lev  ), &
                                      amrex_geom(lev  ), gr_fillPhysicalBC, &
                                      time, &
                                      UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                      amrex_ref_ratio(lev-1), &
                                      amrex_interp_cell_cons, &
                                      lo_bc_amrex, hi_bc_amrex)
    end if

    write(*,'(A,I10,A,I0)') "Created and initialized ", n_blocks, &
                           " blocks on level ", lev + 1

end subroutine gr_initNewLevelCallback

