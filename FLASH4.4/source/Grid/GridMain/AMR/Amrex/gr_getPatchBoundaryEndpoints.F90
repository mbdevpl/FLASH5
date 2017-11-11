!!****if* source/Grid/GridMain/AMR/Amrex/gr_getPatchBoundaryEndpoints
!!
!! NAME
!!
!!  gr_getPatchBoundaryEndpoints
!!
!! SYNOPSIS
!!
!!  call gr_getPatchBoundaryEndpoints()
!!
!! DESCRIPTION 
!!  
!!   volume defined by endPts must contain the volume defined by limits.  They
!!   may be the same volume.
!!
!! ARGUMENTS 
!!  
!!   face - 
!!   axis - 
!!   limits - 0-based, cell-centered, and global
!!   delta - 
!!   endPts - 1-based, cell-centered, and global
!!
!! NOTES
!!
!!
!!  
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_getPatchBoundaryEndpoints(face, axis, limits, delta, endPts)
    use amrex_geometry_module, ONLY : amrex_problo

    use Grid_data,             ONLY : gr_globalDomain
    use Driver_interface,      ONLY : Driver_abortFlash

    implicit none

    integer, intent(IN)  :: face
    integer, intent(IN)  :: axis
    integer, intent(IN)  :: limits(LOW:HIGH, 1:MDIM)
    real,    intent(IN)  :: delta(1:MDIM)
    integer, intent(OUT) :: endPts(LOW:HIGH, 1:MDIM)

    real    :: coord
    integer :: j

    ! Work in 1-based, cell-centered, global indices for FLASH
    
    ! Assume WHOLE_VECTOR along all directions in domain
    endPts(:, :) = 1
    endPts(LOW,  1:NDIM) = limits(LOW,  1:NDIM) + 1
    endPts(HIGH, 1:NDIM) = limits(HIGH, 1:NDIM) + 1

    associate(x0 => amrex_problo(axis), &
              lo => limits(LOW,  axis) + 1, &
              hi => limits(HIGH, axis) + 1, &
              dx => delta(axis))

        endPts(LOW,  axis) = lo
        endPts(HIGH, axis) = lo
        do j = lo, hi
            coord = x0 + (j - 0.5d0)*dx
            if (ABS(coord - gr_globalDomain(face, axis)) < dx) then
                endPts(LOW,  axis) = j - NGUARD + 1
                endPts(HIGH, axis) = j + NGUARD
!                write(*,*) "Found boundary at ", coord, gr_globalDomain(face, axis), dx
                EXIT
            end if
        end do

    end associate

!    write(*,*) "Endpoint Lower = ", endPts(LOW,  :)
!    write(*,*) "Endpoint Upper = ", endPts(HIGH, :)

end subroutine gr_getPatchBoundaryEndpoints

