!!****if* source/Grid/GridMain/AMR/Amrex/gr_untransformBcRegion
!!
!! NAME
!!
!!  gr_untransformBcRegion
!!
!! SYNOPSIS
!!
!!  call gr_untransformBcRegion()
!!
!! DESCRIPTION 
!!  
!!
!!
!! ARGUMENTS 
!!   destData - a pointer to the FAB to which data will be transferred.  Since
!!              it is an AMReX-related element, the spatial indices are 0-based
!!              and have absolute, global indices.  The fourth index is 1-based.
!!   axis - the direction along which the BC were filled (i.e. one of [IJX]AXIS)
!!   limitsGC - the location and size of destSrc given as lower and upper points
!!              in a 1-based cell-centered index space.
!!   region - the data array that contains the data (including the boundary) to
!!            be transferred to destData.  The meaning of the indices is given
!!            in the documentation for XXX.  The indices for this are 1-based 
!!
!!
!! NOTES
!!
!!
!!  
!!***

#include "constants.h"
#include "Flash.h"

! TODO: We want to loop over limits, not endpoints.  We need an offset then for region
subroutine gr_untransformBcRegion(region, axis, boundary, regionSize, destData)
    use amrex_fort_module, ONLY : wp => amrex_real

    implicit none

    real(wp), intent(IN)                         :: region(:, :, :, :)
    integer,  intent(IN)                         :: axis
    integer,  intent(IN)                         :: boundary(LOW:HIGH, 1:MDIM)
    integer,  intent(IN)                         :: regionSize(4)
    real(wp), intent(INOUT), pointer, contiguous :: destData(:, :, :, :)

    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)

    integer :: i, j, k, var
    integer :: n, m

    ! Work in 0-based, cell-centered index space with global indices for AMReX

    lo(:) = 1
    hi(:) = 1
    lo(1:NDIM) = boundary(LOW,  1:NDIM)
    hi(1:NDIM) = boundary(HIGH, 1:NDIM)

    associate (strt   => lo(axis), &
               fin    => hi(axis), &
               nVars  => regionSize(STRUCTSIZE), &
               bcSize => (hi(axis) - lo(axis) + 1))
 
        ! n, m must be 1-based for FLASH / i, j, k must be 0-based for AMReX
        ! var is 1-based for both
        if (axis == IAXIS) then
            do k = lo(KAXIS), hi(KAXIS)
                m = k - lo(KAXIS) + 1

                do j = lo(JAXIS), hi(JAXIS)
                    n = j - lo(JAXIS) + 1
                    do var = 1, nVars
                        destData(strt:fin, j, k, var) = region(1:bcSize, n, m, var)
                    end do
                end do

            end do
        else if (axis == JAXIS) then
            do k = lo(KAXIS), hi(KAXIS)
                m = k - lo(KAXIS) + 1
                
                do i = lo(IAXIS), hi(IAXIS)
                    n = i - lo(IAXIS) + 1
                    do var=1, nVars
                        destData(i, strt:fin, k, var) = region(1:bcSize, n, m, var)
                    end do
                end do

            end do
        else if (axis == KAXIS) then
            do j = lo(JAXIS), hi(JAXIS)
                m = j - lo(JAXIS) + 1
                
                do i = lo(IAXIS), hi(IAXIS)
                    n = i - lo(IAXIS) + 1
                    do var = 1, nVars
                        destData(i, j, strt:fin, var) = region(1:bcSize, n, m, var)
                    end do
                end do
            end do
        end if

    end associate

end subroutine gr_untransformBcRegion

