#include "constants.h"
#include "Flash.h"

subroutine gr_transformBcRegion(srcData, axis, intersection, regionSize, region)
    use amrex_fort_module, ONLY : wp => amrex_real

    implicit none

    real(wp), pointer, contiguous, intent(IN)  :: srcData(:, :, :, :)
    integer,                       intent(IN)  :: axis
    integer,                       intent(IN)  :: intersection(LOW:HIGH, 1:MDIM)
    integer,                       intent(IN)  :: regionSize(4)
    real(wp), pointer, contiguous, intent(OUT) :: region(:, :, :, :)

    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)

    integer :: i, j, k, var
    integer :: n, m

    ! n, m must be 1-based for FLASH / i, j, k must be 0-based for AMReX
    ! var is 1-based for both
    lo(:) = 1
    hi(:) = 1
    lo(1:NDIM) = intersection(LOW,  1:NDIM)
    hi(1:NDIM) = intersection(HIGH, 1:NDIM)

    associate (strt   => lo(axis), &
               fin    => hi(axis), &
               nVars  => regionSize(STRUCTSIZE), &
               bcSize => hi(axis) - lo(axis) + 1)

        ! n, m are 1-based, cell-centered, and local
        if (axis == IAXIS) then
            do k = lo(KAXIS), hi(KAXIS)
                m = k - lo(KAXIS) + 1

                do j = lo(JAXIS), hi(JAXIS)
                    n = j - lo(JAXIS) + 1
                    do var = 1, nVars
                        region(1:bcSize, n, m, var) = srcData(strt:fin, j, k, var)
                    end do
                end do

            end do
        else if (axis == JAXIS) then
            do k = lo(KAXIS), hi(KAXIS)
                m = k - lo(KAXIS) + 1

                do i = lo(IAXIS), hi(IAXIS)
                    n = i - lo(IAXIS) + 1
                    do var = 1, nVars
                        region(1:bcSize, n, m, var) = srcData(i, strt:fin, k, var)
                    end do
                end do

            end do
        else if (axis == KAXIS) then
            do j = lo(JAXIS), hi(JAXIS)
                m = j - lo(JAXIS) + 1

                do i = lo(IAXIS), hi(IAXIS)
                    n = i - lo(IAXIS) + 1
                    do var = 1, nVars
                        region(1:bcSize, n, m, var) = srcData(i, j, strt:fin, var)
                    end do
                end do

            end do
        end if

    end associate

end subroutine gr_transformBcRegion

