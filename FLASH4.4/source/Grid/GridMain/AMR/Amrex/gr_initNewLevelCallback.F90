#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine gr_initNewLevelCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,      ONLY : wp => amrex_real
    use amrex_amr_module,       ONLY : amrex_geom, &
                                       amrex_problo
    use amrex_box_module,       ONLY : amrex_box
    use amrex_boxarray_module,  ONLY : amrex_boxarray
    use amrex_distromap_module, ONLY : amrex_distromap
    use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                       amrex_mfiter_build, &
                                       amrex_mfiter_destroy, &
                                       amrex_multifab_build
    
    use gr_physicalMultifabs,   ONLY : unk, &
                                       facevarx, facevary, facevarz
    use gr_amrexInterface,      ONLY : gr_clearLevelCallback
    use block_metadata,         ONLY : block_metadata_t
    use Simulation_interface,   ONLY : Simulation_initBlock
    use Grid_data,              ONLY : gr_iguard, &
                                       gr_eosModeInit
    use Eos_interface,          ONLY : Eos_wrapped

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
    real(wp)                      :: x = 0.0_wp 
    real(wp)                      :: y = 0.0_wp
    real(wp)                      :: z = 0.0_wp

    integer :: i, j, k

    integer :: n_blocks

#ifdef DEBUG_GRID
    write(*,'(A,A,I2)') "[gr_initNewLevelCallback]", &
                        "            Start Level ", lev + 1
#endif

    ba = pba
    dm = pdm

    call gr_clearLevelCallback(lev)

    ! Create FABS for storing physical data at given level
    call amrex_multifab_build(unk     (lev), ba, dm, NUNK_VARS, gr_iguard)
    ! DEVNOTE: TODO Create there w.r.t. proper face-centered boxes
#if NFACE_VARS > 0
    call amrex_multifab_build(facevarx(lev), ba, dm, NUNK_VARS, gr_iguard)
    call amrex_multifab_build(facevary(lev), ba, dm, NUNK_VARS, gr_iguard)
    call amrex_multifab_build(facevarz(lev), ba, dm, NUNK_VARS, gr_iguard)
#endif

    ! Write initial data across domain at coarsest level
    call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)

    n_blocks = 0
    do while (mfi%next())
        bx = mfi%tilebox()

        ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
        ! Level must be 1-based index and limits/limitsGC must be 1-based also
        ! DEVNOTE: Should we use gr_[ijk]guard here?
        block%level = lev + 1
        block%grid_index = mfi%grid_index()
        block%limits(LOW,  :) = 1
        block%limits(HIGH, :) = 1
        block%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
        block%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
        block%limitsGC(LOW,  :) = 1
        block%limitsGC(HIGH, :) = 1
        block%limitsGC(LOW,  1:NDIM) = block%limits(LOW,  1:NDIM) - NGUARD
        block%limitsGC(HIGH, 1:NDIM) = block%limits(HIGH, 1:NDIM) + NGUARD

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
 
    write(*,'(A,I10,A,I2)') "[gr_initNewLevelCallback]      ", &
                           n_blocks, " new blocks on level", lev + 1

end subroutine gr_initNewLevelCallback

