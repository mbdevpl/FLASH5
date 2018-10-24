!!****if* source/Grid/GridMain/AMR/Amrex/gr_splitFabAtBoundary
!!
!! NAME
!!
!!  gr_splitFabAtBoundary
!!
!! SYNOPSIS
!!
!!  call gr_splitFabAtBoundary(integer(IN)  :: face,
!!                             integer(IN)  :: axis,
!!                             integer(IN)  :: patch(LOW:HIGH, MDIM),
!!                             real(IN)     :: delta(MDIM),
!!                             integer(OUT) :: interior(LOW:HIGH, MDIM),
!!                             integer(OUT) :: guardcells(LOW:HIGH, MDIM),
!!                             logical(OUT) :: hasGC)
!!
!! DESCRIPTION 
!! 
!!   Given a patch and a domain boundary specified by (face, axis), this routine
!!   splits the patch into two regions in the case that the patch straddles the
!!   given boundary.  The routine assumes that the boundary extends the full
!!   length of the patch along the directions parallel to the boundary and the
!!   two regions are those regions on either side of this plane.
!!
!! ARGUMENTS 
!!  
!!   face - specify with a value of LOW or HIGH the boundary of interest
!!   axis - specify with a value of {I,J,K}AXIS the direction to which the
!!          boundary of interest is perpendicular
!!   patch - the specification of the patch through its lower and upper corner
!!           points in a 0-based, cell-centered index space
!!   delta - an array containing the physical cell size
!!   interior - the specification of the region on the side of the boundary that
!!              contains interior cells.  It is defined by its lower and upper 
!!              points in a 0-based, cell-centered index space.
!!   guardcells - the specification of the region on the side of the boundary
!!                that contains guardcells.  It is defined by its lower and
!!                upper points in a 0-based, cell-centered index space.
!!   hasGC - True if the boundary splits the patch into two regions
!!
!! NOTES
!!
!!   If the given patch is located at corner of the domain, the region interior
!!   may contain both interior and guard cells.
!!
!!   If the patch does not straddle the given boundary, then the contents of
!!   interior and guardcells are undefined.
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_splitFabAtBoundary(face, axis, patch, delta, &
                                 interior, guardcells, hasGC)
    use amrex_geometry_module, ONLY : amrex_problo

    use Grid_data,             ONLY : gr_globalDomain

    implicit none

    integer, intent(IN)  :: face
    integer, intent(IN)  :: axis
    integer, intent(IN)  :: patch(LOW:HIGH, 1:MDIM)
    real,    intent(IN)  :: delta(1:MDIM)
    integer, intent(OUT) :: interior(LOW:HIGH, 1:MDIM)
    integer, intent(OUT) :: guardcells(LOW:HIGH, 1:MDIM)
    logical, intent(OUT) :: hasGC

    real    :: coord
    integer :: j

    ! Work in 0-based, cell-centered, global indices for AMReX

    ! Assume boundary extends fully across patch along other directions
    guardcells(:, :) = 1
    guardcells(LOW,  1:NDIM) = patch(LOW,  1:NDIM)
    guardcells(HIGH, 1:NDIM) = patch(HIGH, 1:NDIM)

    interior(:, :) = 1
    interior(LOW,  1:NDIM) = patch(LOW,  1:NDIM)
    interior(HIGH, 1:NDIM) = patch(HIGH, 1:NDIM)

    ! DEV: TODO Improve implementation to allow for irregular shaped
    ! domains defined by Simulation_defineDomain
    associate(x0 => amrex_problo(axis), &
              lo => patch(LOW,  axis), &
              hi => patch(HIGH, axis), &
              dx => delta(axis))

        ! Scan for boundary
        hasGC = .FALSE.
        do j = lo, hi
            coord = x0 + (j + 0.5)*dx
            if (ABS(coord - gr_globalDomain(face, axis)) < dx) then
                if (face == LOW) then
                    guardcells(LOW,  axis) = lo
                    guardcells(HIGH, axis) = j
                    interior(LOW,   axis)  = j + 1
                    interior(HIGH,  axis)  = hi
                else
                    interior(LOW,   axis)  = lo
                    interior(HIGH,  axis)  = j
                    guardcells(LOW,  axis) = j + 1
                    guardcells(HIGH, axis) = hi
                end if

                hasGC = .TRUE.
                RETURN
            end if
        end do
 
    end associate
    
end subroutine gr_splitFabAtBoundary

