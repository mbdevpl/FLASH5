#include "Flash.h"
#include "constants.h"

subroutine sim_writeDataPoints(initData, tileDesc, points, values)
    use flash_tile,     ONLY : flash_tile_t 

    implicit none

    real,                           pointer :: initData(:, :, :, :)
    type(flash_tile_t), intent(IN)          :: tileDesc
    real,               intent(IN)          :: points(:, :)
    real,               intent(IN)          :: values(:)

    real :: deltas(1:MDIM)
    real :: bbox(LOW:HIGH, 1:MDIM)

    real, allocatable :: xcoords(:)
    real, allocatable :: ycoords(:)

    integer :: i, j, p

    allocate(xcoords(tileDesc%limits(LOW, IAXIS):tileDesc%limits(HIGH, IAXIS)), &
             ycoords(tileDesc%limits(LOW, JAXIS):tileDesc%limits(HIGH, JAXIS)))
 
    call tileDesc%deltas(deltas)
    call tileDesc%boundBox(bbox)
    call tileDesc%coordinates(IAXIS, CENTER, .FALSE., xcoords)
    call tileDesc%coordinates(JAXIS, CENTER, .FALSE., ycoords)

    do p = 1, SIZE(points, 1)
        associate(pt   => points(p, :), &
                  r_lo => bbox(LOW,  :), &
                  r_hi => bbox(HIGH, :), &
                  lo   => tileDesc%limits(LOW,  :), &
                  hi   => tileDesc%limits(HIGH, :), &
                  dx   => 0.5 * deltas(IAXIS), &
                  dy   => 0.5 * deltas(JAXIS), &
                  lev  => tileDesc%level)

            if (values(p) == 0.0)   CYCLE

            ! Only continue if point in block
            if (      (r_lo(1) <= pt(1)) .AND. (pt(1) <= r_hi(1)) &
                .AND. (r_lo(2) <= pt(2)) .AND. (pt(2) <= r_hi(2))) then

                ! Search across all cells in box
      cellLoop: do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        if (      (xcoords(i)-dx <= pt(1)) .AND. (pt(1) <= xcoords(i)+dx) &
                            .AND. (ycoords(j)-dy <= pt(2)) .AND. (pt(2) <= ycoords(j)+dy)) then
                            write(*,'(A,I3,A,I3,A)') &
                                  "     Data point contained in cell (", &
                                  i, ",", j, ")"
                            write(*,'(A,F7.5,A,F7.5,A,F7.5)') "        ", &
                                  xcoords(i)-dx, "<=", pt(1), "<=", xcoords(i)+dx
                            write(*,'(A,F7.5,A,F7.5,A,F7.5)') "        ", &
                                  ycoords(j)-dy, "<=", pt(2), "<=", ycoords(j)+dy

                            ! Using 1-based FLASH level indexing
                            !
                            ! Given values are integers that indicate level of
                            ! refinement to achieve at point.  However, this
                            ! simulation averages data across levels and
                            ! therefore we must store densities.
                            initData(i, j, 1, :) = values(p) * 4.0**(lev-4)

                            exit cellLoop
                        end if
                    end do
                end do cellLoop
            
            end if
        
        end associate
    end do

    deallocate(xcoords, ycoords)

end subroutine sim_writeDataPoints

