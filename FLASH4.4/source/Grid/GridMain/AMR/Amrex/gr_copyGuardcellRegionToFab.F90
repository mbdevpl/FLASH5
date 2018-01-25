!!****if* source/Grid/GridMain/AMR/Amrex/gr_copyGuardcellRegionToFab
!!
!! NAME
!!
!!  gr_copyGuardcellRegionToFab
!!
!! SYNOPSIS
!!
!!  call gr_copyGuardcellRegionToFab(real(IN)   :: region(:,:,:,:),
!!                                  integer(IN) :: face,
!!                                  integer(IN) :: axis,
!!                                  integer(IN) :: guardcells(LOW:HIGH, MDIM),
!!                                  integer(IN) :: scomp,
!!                                  integer(IN) :: ncomp,
!!                                  real(INOUT) :: fab(:,:,:,:))
!!
!! DESCRIPTION 
!!  
!!  This routine is used to copy guardcell data from a special data structure
!!  filled by Grid_bcApplyToRegion into the AMReX FAB that manages physical
!!  data.
!!
!! ARGUMENTS
!!
!!  region - Data structure that sources GC data.  Please see 
!!           Grid_bcApplyToRegion.
!!  face - specify with a value of LOW or HIGH the boundary of interest
!!  axis - specify with a value of {I,J,K}AXIS the direction to which the
!!         boundary of interest is perpendicular
!!  guardcells - the specification of the GC region into which data shall be
!!               copied.  It is defined by its lower and upper points in a
!!               0-based, cell-centered index space.
!!  scomp - the 1-based index of the first physical quantity to copy
!!  ncomp - the number of physical quantities to copy starting from scomp
!!  fab - a pointer to the FAB to which data will be transferred.  Since
!!        it is an AMReX-related element, the spatial indices are 0-based
!!        and have absolute, global indices.  The fourth index is 1-based.
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_copyGuardcellRegionToFab(region, face, axis, guardcells, &
                                       scomp, ncomp, fab)
    use amrex_fort_module, ONLY : wp => amrex_real

#ifdef DEBUG_GRID
    use Driver_interface,  ONLY : Driver_abortFlash
#endif

    implicit none

    real(wp), intent(IN),    pointer, contiguous :: region(:, :, :, :)
    integer,  intent(IN)                         :: face
    integer,  intent(IN)                         :: axis
    integer,  intent(IN)                         :: guardcells(LOW:HIGH, 1:MDIM)
    integer,  intent(IN)                         :: scomp
    integer,  intent(IN)                         :: ncomp
    real(wp), intent(INOUT), pointer, contiguous :: fab(:, :, :, :)

    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)

    integer :: i, j, k, var
    integer :: n, m

    integer :: rStrt, rFin

    ! n, m must be 1-based, cell-centered for FLASH and local for region
    ! i, j, k must be 0-based, cell-centered for AMReX and global for FAB
    ! var is 1-based for both

    ! Assume boundary extends fully across patch along other directions
    lo(:) = 1
    hi(:) = 1
    lo(1:NDIM) = guardcells(LOW,  1:NDIM)
    hi(1:NDIM) = guardcells(HIGH, 1:NDIM)
    
    associate (strt  => lo(axis), &
               fin   => hi(axis), &
               width => (hi(axis) - lo(axis) + 1))
   
#ifdef DEBUG_GRID
        if (width > NGUARD) then
            call Driver_abortFlash("[gr_copyGuardcellsToFab] Given patch is too wide")
        end if
#endif

        ! Determine offset to GC along BC direction
        ! NOTE: The width of the intersection of patch and GC region
        ! can be can be less than NGUARD
        if (face == LOW) then
            rFin  = NGUARD
            rStrt = rFin  - (width - 1)
        else
            rStrt = NGUARD + 1
            rFin  = rStrt + (width - 1)
        end if
 
        if (axis == IAXIS) then
            do k = lo(KAXIS), hi(KAXIS)
                m = k - lo(KAXIS) + 1

                do j = lo(JAXIS), hi(JAXIS)
                    n = j - lo(JAXIS) + 1
                    do var = 1, ncomp
                        fab(strt:fin, j, k, var+scomp-1) = region(rStrt:rFin, n, m, var)
                    end do
                end do

            end do
        else if (axis == JAXIS) then
            do k = lo(KAXIS), hi(KAXIS)
                m = k - lo(KAXIS) + 1
                
                do i = lo(IAXIS), hi(IAXIS)
                    n = i - lo(IAXIS) + 1
                    do var = 1, ncomp
                        fab(i, strt:fin, k, var+scomp-1) = region(rStrt:rFin, n, m, var)
                    end do
                end do

            end do
        else if (axis == KAXIS) then
            do j = lo(JAXIS), hi(JAXIS)
                m = j - lo(JAXIS) + 1
                
                do i = lo(IAXIS), hi(IAXIS)
                    n = i - lo(IAXIS) + 1
                    do var = 1, ncomp
                        fab(i, j, strt:fin, var+scomp-1) = region(rStrt:rFin, n, m, var)
                    end do
                end do
            end do
        end if

    end associate

end subroutine gr_copyGuardcellRegionToFab

