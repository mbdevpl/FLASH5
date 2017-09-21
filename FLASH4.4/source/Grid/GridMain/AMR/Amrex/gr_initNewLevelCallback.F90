subroutine gr_initNewLevelCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,      ONLY : wp => amrex_real
    use amrex_amr_module,       ONLY : amrex_geom, &
                                       amrex_problo
    use amrex_box_module,       ONLY : amrex_box
    use amrex_boxarray_module,  ONLY : amrex_boxarray
    use amrex_distromap_module, ONLY : amrex_distromap
    use amrex_parallel_module,  ONLY : amrex_parallel_myproc
    use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                       amrex_mfiter_build, &
                                       amrex_mfiter_destroy, &
                                       amrex_multifab_build
    
    use gr_physicalMultifabs,   ONLY : unk, &
                                       facevarx, facevary, facevarz
    use amrex_interfaces,       ONLY : gr_clearLevelCallback
    use block_metadata,         ONLY : block_metadata_t
    use Simulation_interface,   ONLY : Simulation_initBlock
    use Grid_data,              ONLY : gr_iguard

    implicit none

#include "constants.h"
#include "Flash.h"

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

    integer :: rank = 0

    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
 
    rank = amrex_parallel_myproc()
    write(*,*) "[gr_initNewLevelCallback] Start Level ", lev + 1
 
    ba = pba
    dm = pdm

    call gr_clearLevelCallback(lev)

    ! Create FABS for storing physical data at coarsest level
    call amrex_multifab_build(unk     (lev), ba, dm, NUNK_VARS, gr_iguard)
    ! DEVNOTE: TODO Create test wrt proper face-centered boxes
    call amrex_multifab_build(facevarx(lev), ba, dm, NUNK_VARS, gr_iguard)
    call amrex_multifab_build(facevary(lev), ba, dm, NUNK_VARS, gr_iguard)
    call amrex_multifab_build(facevarz(lev), ba, dm, NUNK_VARS, gr_iguard)

    ! Write initial data across domain at coarsest level
    call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)

    do while (mfi%next())
        bx = mfi%tilebox()

        ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
        ! Level must be 1-based index and limits/limitsGC must be 1-based also
        ! DEVNOTE: Should we use gr_[ijk]guard here?
        block%level = lev + 1
        ! DEVNOTE: TODO Get grid_index from mfi
        block%grid_index = -1
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
        nullify(initData)
    end do

    call amrex_mfiter_destroy(mfi)
 
    write(*,*) "[gr_initNewLevelCallback] Finished Level ", lev + 1
end subroutine gr_initNewLevelCallback

