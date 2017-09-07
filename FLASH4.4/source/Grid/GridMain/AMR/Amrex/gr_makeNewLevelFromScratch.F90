subroutine gr_makeNewLevelFromScratch(lev, time, pba, pdm) bind(c)
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
    use amrex_interfaces,       ONLY : gr_clearLevel
    use block_metadata,         ONLY : block_metadata_t

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
 
    integer :: lev_flash = -1

    ! C++ AMReX uses zero-based level index set, but FLASH uses 1-based set
    lev_flash = lev + 1 

    rank = amrex_parallel_myproc()
    write(*,*) "[Rank ", rank, "] gr_makeNewLevelFromScratch - Start Level ", &
                lev_flash
 
    ba = pba
    dm = pdm

    call gr_clearLevel(lev)

    ! Create FABS for storing physical data at coarsest level
    call amrex_multifab_build(unk(lev_flash),      ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevarx(lev_flash), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevary(lev_flash), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevarz(lev_flash), ba, dm, NUNK_VARS, NGUARD)

    ! Write initial data across domain at coarsest level
    call amrex_mfiter_build(mfi, unk(lev_flash))

    do while (mfi%next())
        bx = mfi%tilebox()
        initData => unk(lev_flash)%dataptr(mfi)

        ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
        block%level = lev
        block%grid_index = -1
        block%limits(LOW,  :) = bx%lo
        block%limits(HIGH, :) = bx%hi
        block%limitsGC(LOW, :) = bx%lo - NGUARD
        block%limitsGC(HIGH, :) = bx%hi + NGUARD

        !  We need to zero data in case we reuse blocks from previous levels
        !  but don't initialize all data in Simulation_initBlock... in particular
        !  the total vs. internal energies can cause problems in the eos call that 
        !  follows.
        initData = 0.0
        call Simulation_initBlock(initData, block)
        nullify(initData)
    end do

    call amrex_mfiter_destroy(mfi)
 
    write(*,*) "[Rank ", rank, "] gr_makeNewLevelFromScratch - End Level ", &
                lev_flash
end subroutine gr_makeNewLevelFromScratch 

