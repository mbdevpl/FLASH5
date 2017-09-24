subroutine sim_writeDataPoints(initData, blockDesc, points, values)
    use block_metadata, ONLY : block_metadata_t
    use Grid_interface, ONLY : Grid_getDeltas, &
                               Grid_getBlkBoundBox, &
                               Grid_getSingleCellCoords

    implicit none

    real,                   intent(IN), pointer :: initData(:, :, :, :)
    type(block_metadata_t), intent(IN)          :: blockDesc
    real,                   intent(IN)          :: points(:, :)
    real,                   intent(IN)          :: values(:)

#include "Flash.h"
#include "constants.h"

    real    :: deltas(1:MDIM)
    real    :: bbox(LOW:HIGH, 1:MDIM)
    real    :: cd(1:MDIM)
    integer :: idx(1:MDIM)

    integer :: i, j, p

    deltas = 0.0d0
    call Grid_getDeltas(blockDesc%level, deltas)
   
    bbox(:, :) = 0.0d0
    call Grid_getBlkBoundBox(blockDesc, bbox)

    do p = 1, SIZE(points, 1)
        associate(pt   => points(p, :), &
                  r_lo => bbox(LOW,  :), &
                  r_hi => bbox(HIGH, :), &
                  lo   => blockDesc%limits(LOW,  :), &
                  hi   => blockDesc%limits(HIGH, :), &
                  dx   => 0.5d0 * deltas(IAXIS), &
                  dy   => 0.5d0 * deltas(JAXIS), &
                  lev  => blockDesc%level)

            if (values(p) == 0.0d0)   CYCLE
            
            ! Only continue if point in block
            if (      (r_lo(1) <= pt(1)) .AND. (pt(1) <= r_hi(1)) &
                .AND. (r_lo(2) <= pt(2)) .AND. (pt(2) <= r_hi(2))) then
                
                ! Search across all cells in box
      cellLoop: do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        ! Use local index to get center coordinate of cell
                        idx = [i-lo(IAXIS)+1, j-lo(JAXIS)+1, 1]
                        call Grid_getSingleCellCoords(idx, blockDesc, CENTER, &
                                                      INTERIOR, cd)

                        if (      (cd(1)-dx <= pt(1)) .AND. (pt(1) <= cd(1)+dx) &
                            .AND. (cd(2)-dy <= pt(2)) .AND. (pt(2) <= cd(2)+dy)) then
                            write(*,'(A,I2,A,I2,A)') &
                                  "     Data point contained in cell (", &
                                  i, ",", j, ")"
                            write(*,'(A,F7.5,A,F7.5,A,F7.5)') "        ", &
                                  cd(1)-dx, "<=", pt(1), "<=", cd(1)+dx
                            write(*,'(A,F7.5,A,F7.5,A,F7.5)') "        ", &
                                  cd(2)-dy, "<=", pt(2), "<=", cd(2)+dy

                            ! Using 1-based FLASH level indexing
                            !
                            ! Given values are integers that indicate level of
                            ! refinement to achieve at point.  However, this
                            ! simulation averages data across levels and
                            ! therefore we must store densities.
                            initData(i, j, 1, :) = values(p) * 4.0d0**(lev-4)

                            exit cellLoop
                        end if
                    end do
                end do cellLoop
            
            end if
        
        end associate
    end do

end subroutine sim_writeDataPoints

