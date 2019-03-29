!!****if* source/Grid/GridMain/gr_getRegionDataCoordinates.F90
!!
!! NAME
!!  gr_getRegionDataCoordinates
!!
!! SYNOPSIS
!!  call Grid_bcApplyToRegionSpecialized(integer(IN)  :: level, 
!!                                       integer(IN)  :: gridDataStruct,
!!                                       integer(IN)  :: axis,
!!                                       integer(IN)  :: axis2,
!!                                       integer(IN)  :: axis3,
!!                                       integer(IN)  :: regionSize(REGION_DIM),
!!                                       integer(IN)  :: endPoints(LOW:HIGH, MDIM),
!!                                       real(OUT)    :: coordinates(regionSize(BC_DIR),
!!                                                                   regionSize(SECOND_DIR), 
!!                                                                   regionSize(THIRD_DIR),
!!                                                                   1:MDIM)) 
!!
!! DESCRIPTION
!!  This is a helper routine for those versions of Grid_bcApplyToRegionSpecialized
!!  that require knowledge of the physical coordinates of all mesh elements in
!!  the given regionData on which they are to operate.
!!
!! ARGUMENTS
!!  See ARGUMENTS section of documentation for Grid_bcApplyToRegion.
!!
!!***

#include "constants.h"

subroutine gr_getRegionDataCoordinates(level, gridDataStruct, &
                                       axis, axis2, axis3, &
                                       regionSize, endPoints, &
                                       coordinates)
  use Grid_data,        ONLY : gr_globalDomain 
  use Grid_interface,   ONLY : Grid_getDeltas

  implicit none

  integer,                intent(IN)  :: level
  integer,                intent(IN)  :: gridDataStruct
  integer,                intent(IN)  :: axis
  integer,                intent(IN)  :: axis2
  integer,                intent(IN)  :: axis3
  integer,                intent(IN)  :: regionSize(REGION_DIM)
  integer,                intent(IN)  :: endPoints(LOW:HIGH, 1:MDIM)
  real,                   intent(OUT) :: coordinates(regionSize(BC_DIR),     &
                                                     regionSize(SECOND_DIR), &
                                                     regionSize(THIRD_DIR),  &
                                                     1:MDIM)

  integer :: i, j, k
  integer :: i_global, j_global, k_global
  real    :: i_shift, j_shift, k_shift

  real    :: x0(1:MDIM)
  real    :: dx(1:MDIM)

  ! Assume cell-centered index space
  i_shift = 0.5
  j_shift = 0.5
  k_shift = 0.5

  ! Adjust for face-centered index space if necessary.
  ! Assume no edge-centered or nodal.
  if      (     ((gridDataStruct == FACEX) .AND. (axis == IAXIS)) &
           .OR. ((gridDataStruct == FACEY) .AND. (axis == JAXIS)) &
           .OR. ((gridDataStruct == FACEZ) .AND. (axis == KAXIS)) ) then
    i_shift = 1.0
  else if (     ((gridDataStruct == FACEX) .AND. (axis2 == IAXIS)) &
           .OR. ((gridDataStruct == FACEY) .AND. (axis2 == JAXIS)) &
           .OR. ((gridDataStruct == FACEZ) .AND. (axis2 == KAXIS)) ) then
    j_shift = 1.0
  else if (     ((gridDataStruct == FACEX) .AND. (axis3 == IAXIS)) &
           .OR. ((gridDataStruct == FACEY) .AND. (axis3 == JAXIS)) &
           .OR. ((gridDataStruct == FACEZ) .AND. (axis3 == KAXIS)) ) then
    k_shift = 1.0
  end if

  x0(:) = gr_globalDomain(LOW, :)
  call Grid_getDeltas(level, dx)
  do     k = 1, regionSize(THIRD_DIR)
    do   j = 1, regionSize(SECOND_DIR)
      do i = 1, regionSize(BC_DIR)
        ! Determine location of current cell in global index space
        i_global = endPoints(LOW, axis ) + i - 1
        j_global = endPoints(LOW, axis2) + j - 1
        k_global = endPoints(LOW, axis3) + k - 1

        ! Adjust to mesh element of index space and to physical coordinates
        coordinates(i, j, k, axis ) = x0(axis ) + (i_global - i_shift)*dx(axis )
        coordinates(i, j, k, axis2) = x0(axis2) + (j_global - j_shift)*dx(axis2)
        coordinates(i, j, k, axis3) = x0(axis3) + (k_global - k_shift)*dx(axis3)
      end do
    end do
  end do
end subroutine gr_getRegionDataCoordinates

