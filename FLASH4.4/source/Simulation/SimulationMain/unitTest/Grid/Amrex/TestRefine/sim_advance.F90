subroutine sim_advance(step, points, values, set_msg, leaf_msg, block_count)
    use Grid_interface,       ONLY : Grid_updateRefinement, &
                                     Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                     Grid_averageDownLevels
    use Driver_interface,     ONLY : Driver_abortFlash
    use block_iterator,       ONLY : block_iterator_t 
    use block_metadata,       ONLY : block_metadata_t 
    use amrex_interfaces,     ONLY : gr_getFinestLevel
    use sim_interface,        ONLY : sim_writeDataPoints, &
                                     sim_printLeaves

    implicit none

#include "constants.h"
    
    integer,      intent(IN)    :: step
    real,         intent(IN)    :: points(:, :)
    real,         intent(IN)    :: values(:)
    character(*), intent(IN)    :: set_msg
    character(*), intent(IN)    :: leaf_msg
    integer,      intent(INOUT) :: block_count(4)

    real, contiguous, pointer :: solnData(:,:,:,:)
   
    type(block_iterator_t) :: itor
    type(block_metadata_t) :: blockDesc

    logical :: gridChanged
    integer :: finest_level
    integer :: lev

    !!!!! ZERO DATA AND WRITE GIVEN POINTS
    write(*,*)
    write(*,'(A,I2,A,I2)') "ADVANCE STEPS ", step, "/", step+1
    write(*,'(A)') "--------------------------------------------------------------"
    write(*,*) set_msg
    write(*,*)

    ! Write to leaf blocks first.  AMReX level indexing is 0-based
    call gr_getFinestLevel(finest_level)
    do lev = 1, finest_level
        itor = block_iterator_t(LEAF, level=lev)
        do while (itor%is_valid())
            call itor%blkMetaData(blockDesc)
            call Grid_getBlkPtr(blockDesc, solnData)

            solnData = 0.0d0
            call sim_writeDataPoints(solnData, blockDesc,  points, values)

            call Grid_releaseBlkPtr(blockDesc, solnData)
            call itor%next()
        end do
    end do

    ! Propogate leaf data to coarse, non-leaf blocks
    call Grid_averageDownLevels

    !!!!! REFINE/DEREFINE BASED ON NEW DATA
    write(*,*)
    write(*,*) "EXECUTING REFINE/DEREFINE"
    write(*,*)

    ! Should only refine on every other step (nrefs = 2)
    gridChanged = .FALSE.
    call Grid_updateRefinement(step,   DBLE(step  ), gridChanged) 
    if (gridChanged) then
        call Driver_abortFlash("[sim_advance] Should not refine on odd steps")
    end if

    call Grid_updateRefinement(step+1, DBLE(step+1), gridChanged) 
    if (.NOT. gridChanged) then
        call Driver_abortFlash("[sim_advance] Should refine on even steps")
    end if

    write(*,*)
    call sim_printLeaves(leaf_msg, block_count)
end subroutine sim_advance

