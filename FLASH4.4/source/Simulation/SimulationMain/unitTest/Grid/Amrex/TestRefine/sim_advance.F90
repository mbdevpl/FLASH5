#include "constants.h"
 
subroutine sim_advance(step, points, values, set_msg, leaf_msg)
    use Grid_interface,       ONLY : Grid_updateRefinement, &
                                     Grid_getTileIterator, &
                                     Grid_releaseTileIterator
    use gr_amrexInterface,    ONLY : gr_restrictAllLevels
    use flash_iterator,       ONLY : flash_iterator_t 
    use flash_tile,           ONLY : flash_tile_t 
    use Driver_interface,     ONLY : Driver_abortFlash
    use sim_interface,        ONLY : sim_writeDataPoints, &
                                     sim_collectLeaves, &
                                     sim_printLeaves

    implicit none

    integer,      intent(IN) :: step
    real,         intent(IN) :: points(:, :)
    real,         intent(IN) :: values(:)
    character(*), intent(IN) :: set_msg
    character(*), intent(IN) :: leaf_msg

    real, contiguous, pointer :: solnData(:,:,:,:)
 
    type(flash_iterator_t) :: itor
    type(flash_tile_t)     :: tileDesc

    logical :: gridChanged

    nullify(solnData)

    !!!!! ZERO DATA AND WRITE GIVEN POINTS
    write(*,*)
    write(*,'(A,I2,A,I2)') "ADVANCE STEPS ", step, "/", step+1
    write(*,'(A)') "--------------------------------------------------------------"
    write(*,*) set_msg
    write(*,*)

    ! Write to leaf blocks first.  AMReX level indexing is 0-based
    call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
    do while (itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnData, CENTER)

        solnData(:,:,:,:) = 0.0
        call sim_writeDataPoints(solnData, tileDesc, points, values)

        call tileDesc%releaseDataPtr(solnData, CENTER)
        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    ! Propogate leaf data to coarse, non-leaf blocks
    call gr_restrictAllLevels(CENTER, convertPtoC=.FALSE., convertCtoP=.FALSE.)

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

    call sim_collectLeaves
end subroutine sim_advance

