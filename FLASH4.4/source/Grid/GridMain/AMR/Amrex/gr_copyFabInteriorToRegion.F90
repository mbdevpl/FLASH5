!!****if* source/Grid/GridMain/AMR/Amrex/gr_copyFabInteriorToRegion
!!
!! NAME
!!  gr_copyFabInteriorToRegion
!!
!! SYNOPSIS
!!  call gr_copyFabInteriorToRegion(real(IN)    :: fab(:,:,:,:),
!!                                  integer(IN) :: gds,
!!                                  integer(IN) :: face,
!!                                  integer(IN) :: axis,
!!                                  integer(IN) :: interior(LOW:HIGH, MDIM),
!!                                  integer(IN) :: scomp,
!!                                  integer(IN) :: ncomp,
!!                                  real(INOUT) :: region(:,:,:,:))
!!
!! DESCRIPTION 
!!  This routine is used to populate a special data structure with interior data
!!  for the purpose of applying BCs to the structure with the 
!!  Grid_bcApplyToRegion routine.  It copies the interior data (see 
!!  gr_splitFabAtBoundary) from the source fab into this special data structure.
!!  If the given index space is face-centered and the fab contains face data
!!  defined on the boundary, then this data is also copied into the special
!!  data structure.
!!
!!  This routine assumes that the fab, interior, and region data structures are
!!  all specified with respect to the same index space.
!!
!! ARGUMENTS 
!!  fab -  an AMReX FAB containing the interior data.  The spatial indices are
!!         with respect to a 0-based index space.  The index
!!         over physical quantities is 1-based.
!!  gds  - the grid data structure associated with fab, interior, and region.
!!         Acceptable values are CENTER and FACE[XYZ].
!!  face - specify with a value of LOW or HIGH the boundary of interest.
!!  axis - specify with a value of [IJK]AXIS the direction to which the
!!         boundary of interest is perpendicular.
!!  interior - the specification of the region from which data shall be copied. 
!!             It is defined by its lower and upper points in a 0-based
!!             index space.
!!  scomp - the 1-based index of the first physical quantity to copy
!!  ncomp - the number of physical quantities to copy starting from scomp
!!  region - Data structure that receives data.  Please see
!!           Grid_bcApplyToRegion.
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_copyFabInteriorToRegion(fab, gds, face, axis, interior, &
                                      scomp, ncomp, region)
    use amrex_fort_module, ONLY : wp => amrex_real

    use Driver_interface,  ONLY : Driver_abortFlash

    implicit none

    real(wp), pointer, contiguous, intent(IN)    :: fab(:, :, :, :)
    integer,                       intent(IN)    :: gds
    integer,                       intent(IN)    :: face
    integer,                       intent(IN)    :: axis
    integer,                       intent(IN)    :: interior(LOW:HIGH, 1:MDIM)
    integer,                       intent(IN)    :: scomp
    integer,                       intent(IN)    :: ncomp
    real(wp), pointer, contiguous, intent(INOUT) :: region(:, :, :, :)

    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)

    integer :: i, j, k, var
    integer :: n, m

    integer :: rStrt, rFin
    integer :: offset
    integer :: rOffset

    ! n, m must be 1-based for FLASH and local for region
    ! i, j, k must be 0-based for AMReX and global for FAB
    ! var is 1-based for both

    if ((gds /= CENTER) .AND. &
        (gds /= FACEX)  .AND. (gds /= FACEY) .AND. (gds /= FACEZ)) then
        call Driver_abortFlash("[gr_copyFabInteriorToRegion] " // &
                               "GDS must be cell- or face-centered")
    end if

    ! Assume boundary extends fully across patch along other directions
    lo(:) = 1
    hi(:) = 1
    lo(1:NDIM) = interior(LOW,  1:NDIM)
    hi(1:NDIM) = interior(HIGH, 1:NDIM)

    ! Assume that we have cell centers along the BC axis.
    ! Else, we have face centers and must grow by one.
    offset = NGUARD - 1
    rOffset = 0
    if (     ((gds == FACEX) .AND. (axis == IAXIS)) &
        .OR. ((gds == FACEY) .AND. (axis == JAXIS)) &
        .OR. ((gds == FACEZ) .AND. (axis == KAXIS))) then
        offset = NGUARD
        rOffset = 1
    end if

    ! Only need from FAB interior closest NGUARD cells along BC direction
    if (face == LOW) then
        lo(axis) = interior(LOW, axis)
        hi(axis) = interior(LOW, axis)  + offset
        ! Skip over guardcells
        rStrt = NGUARD + 1
        rFin  = 2*NGUARD + rOffset
    else 
        lo(axis) = interior(HIGH, axis) - offset
        hi(axis) = interior(HIGH, axis)
        ! Ignore guardcells at end
        rStrt = 1
        rFin  = NGUARD + rOffset
    end if

    associate (strt   => lo(axis), &
               fin    => hi(axis))

        if (axis == IAXIS) then
            do        var = 1, ncomp
                do      k = lo(KAXIS), hi(KAXIS)
                        m = k - lo(KAXIS) + 1
                    do  j = lo(JAXIS), hi(JAXIS)
                        n = j - lo(JAXIS) + 1
                        region(rStrt:rFin, n, m, var) = fab(strt:fin, j, k, var+scomp-1)
                    end do
                end do
            end do
        else if (axis == JAXIS) then
            do        var = 1, ncomp
                do      k = lo(KAXIS), hi(KAXIS)
                        m = k - lo(KAXIS) + 1
                    do  i = lo(IAXIS), hi(IAXIS)
                        n = i - lo(IAXIS) + 1
                        region(rStrt:rFin, n, m, var) = fab(i, strt:fin, k, var+scomp-1)
                    end do
                end do
            end do
        else if (axis == KAXIS) then
            do        var = 1, ncomp
                do      j = lo(JAXIS), hi(JAXIS)
                        m = j - lo(JAXIS) + 1
                    do  i = lo(IAXIS), hi(IAXIS)
                        n = i - lo(IAXIS) + 1
                        region(rStrt:rFin, n, m, var) = fab(i, j, strt:fin, var+scomp-1)
                    end do
                end do
            end do
        end if

    end associate

end subroutine gr_copyFabInteriorToRegion

