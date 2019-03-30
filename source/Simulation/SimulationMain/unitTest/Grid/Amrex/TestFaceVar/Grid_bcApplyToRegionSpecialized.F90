!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFaceVar/Grid_bcApplyToRegionSpecialized
!!
!! NAME
!!  Grid_bcApplyToRegionSpecialized
!!
!! SYNOPSIS
!!  call Grid_bcApplyToRegionSpecialized(integer(IN)  :: bcType, 
!!                                       integer(IN)  :: gridDataStruct,
!!                                       integer(IN)  :: level,
!!                                       integer(IN)  :: guard,
!!                                       integer(IN)  :: axis,
!!                                       integer(IN)  :: face,
!!                                       real(INOUT)  :: regionData(regionSize(BC_DIR),
!!                                                                  regionSize(SECOND_DIR), 
!!                                                                  regionSize(THIRD_DIR),
!!                                                                  regionSize(STRUCTSIZE))
!!                                       integer(IN)  :: regionSize(REGION_DIM),
!!                                       logical(IN)  :: mask(STRUCTSIZE),
!!                                       logical(OUT) :: applied,
!!                                       integer(IN)  :: secondDir,
!!                                       integer(IN)  :: thirdDir,
!!                                       integer(IN)  :: endPoints(LOW:HIGH, MDIM),
!!                             optional, integer(IN)  :: idest)
!!
!! DESCRIPTION
!!  This is a custom routine for assigning data to the face-centered data that
!!  is defined outside of the physical domain and its boundary.  As this routine
!!  is called by gr_fillPhysicalBC, we can therefore test correct functionality
!!  of different aspects of this routine.
!!
!!  All data values defined in the interior of the physical domain are
!!  untouched.  All those face-centered elements on the domain boundary
!!  have their values set to zero.  All face-centered values defined outside
!!  the domain and its boudary are set to 
!!   - 1.1 for the first face variable and
!!   - 2.2 for the second face variable.
!!
!! ARGUMENTS
!!  See ARGUMENTS section of documentation for Grid_bcApplyToRegion.
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_bcApplyToRegionSpecialized(bcType, gridDataStruct, level, &
                                           guard, axis, face, &
                                           regionData, regionSize, &
                                           mask, applied, &
                                           axis2, axis3, &
                                           endPoints, idest)
  use Grid_data,        ONLY : gr_globalDomain 
  use Grid_interface,   ONLY : Grid_getDeltas
  use gr_interface,     ONLY : gr_getRegionDataCoordinates

  implicit none

  integer,                intent(IN)           :: bcType
  integer,                intent(IN)           :: gridDataStruct
  integer,                intent(IN)           :: level
  integer,                intent(IN)           :: guard
  integer,                intent(IN)           :: axis
  integer,                intent(IN)           :: face
  integer,                intent(IN)           :: regionSize(REGION_DIM)
  real,                   intent(INOUT)        :: regionData(regionSize(BC_DIR),     &
                                                             regionSize(SECOND_DIR), &
                                                             regionSize(THIRD_DIR),  &
                                                             regionSize(STRUCTSIZE))
  logical,                intent(IN)           :: mask(regionSize(STRUCTSIZE))
  logical,                intent(OUT)          :: applied
  integer,                intent(IN)           :: axis2, axis3
  integer,                intent(IN)           :: endPoints(LOW:HIGH, MDIM)
  integer,                intent(IN), OPTIONAL :: idest

  integer :: i, j, k, var
  logical :: isOnBoundary
  logical :: containsBoundary(LOW:HIGH, 1:MDIM)

  real    :: x0(1:MDIM)
  real    :: blkBounds(LOW:HIGH, 1:MDIM)
  real    :: deltas(1:MDIM)
  real, allocatable :: coordinates(:,:,:,:)

#ifdef DEBUG_GRID
  write(*,*) "Specialized BC Handling"
  write(*,*) "---------------------------------------------------"
  if      (axis == IAXIS) then
    write(*,*) "Axis   = x-axis"
  else if (axis == JAXIS) then
    write(*,*) "Axis   = y-axis"
  else
    write(*,*) "Axis   = z-axis"
  end if

  if (face == LOW) then
    write(*,*) "Face   = LOW"
  else
    write(*,*) "Face   = HIGH"
  end if

  if      (gridDataStruct == CENTER) then
    write(*,*) "Index Space = CENTER"
  else if (gridDataStruct == FACEX) then
    write(*,*) "Index Space = FACEX"
  else if (gridDataStruct == FACEY) then
    write(*,*) "Index Space = FACEY"
  else if (gridDataStruct == FACEZ) then
    write(*,*) "Index Space = FACEZ"
  end if
  write(*,*) "NGUARD = ", guard

  write(*,*) "endPoints Low:       ", endPoints(LOW, :)
  write(*,*) "endPoints High:      ", endPoints(HIGH, :)
#endif

  call Grid_getDeltas(level, deltas)

  blkBounds(:, :) = 1.0
  x0(:) = gr_globalDomain(LOW, :)
  associate(lo   => endPoints(LOW,  :), &
            hi   => endPoints(HIGH, :))
    blkBounds(LOW,  1:NDIM) = x0(1:NDIM) + (lo(1:NDIM) - 1)*deltas(1:NDIM)
    blkBounds(HIGH, 1:NDIM) = x0(1:NDIM) + (hi(1:NDIM)    )*deltas(1:NDIM)
  end associate 

  associate(lo => blkBounds(LOW,  :), &
            hi => blkBounds(HIGH, :))
    ! DEV: TODO These should be made inexact
    containsBoundary(:, :) = .FALSE.
    containsBoundary(LOW,  IAXIS) =       (lo(IAXIS) <= gr_globalDomain(LOW,  IAXIS)) &
                                    .AND. (gr_globalDomain(LOW,  IAXIS) <= hi(IAXIS))
    containsBoundary(HIGH, IAXIS) =       (lo(IAXIS) <= gr_globalDomain(HIGH, IAXIS)) &
                                    .AND. (gr_globalDomain(HIGH, IAXIS) <= hi(IAXIS))
    containsBoundary(LOW,  JAXIS) =       (lo(JAXIS) <= gr_globalDomain(LOW,  JAXIS)) &
                                    .AND. (gr_globalDomain(LOW,  JAXIS) <= hi(JAXIS))
    containsBoundary(HIGH, JAXIS) =       (lo(JAXIS) <= gr_globalDomain(HIGH, JAXIS)) &
                                    .AND. (gr_globalDomain(HIGH, JAXIS) <= hi(JAXIS))
  end associate

  isOnBoundary =       (((gridDataStruct == FACEX) .AND. (axis == IAXIS)) &
                 .OR.   ((gridDataStruct == FACEY) .AND. (axis == JAXIS)) &
                 .OR.   ((gridDataStruct == FACEZ) .AND. (axis == KAXIS))) &
                 .AND. containsBoundary(face, axis)

  ! Update data on boundary to see if it is being correctly transferred
  ! back to data structures
  if (isOnBoundary) then

    ! Set all face-centered values on domain boundary to zero.  Those
    ! coplanar to, but not on the boundary are considered guardcells
    if (face == LOW) then
      i = guard+1
    else
      i = regionSize(BC_DIR) - NGUARD
    end if

    allocate(coordinates(SIZE(regionData, 1), &
                         SIZE(regionData, 2), &
                         SIZE(regionData, 3), &
                         1:MDIM))

    call gr_getRegionDataCoordinates(level, gridDataStruct, &
                                     axis, axis2, axis3, &
                                     regionSize, endPoints, &
                                     coordinates)

    do     var = 1, regionSize(STRUCTSIZE)
      do     k = 1, regionSize(THIRD_DIR)
        do   j = 1, regionSize(SECOND_DIR)
          ! By construction, the face is coplanar to boundary.  Check if
          ! it is actually on boundary.
          if (      (gr_globalDomain(LOW,  axis2) <= coordinates(i,j,k,axis2)) &
              .AND. (coordinates(i,j,k,axis2) <= gr_globalDomain(HIGH, axis2)) &
              .AND. (gr_globalDomain(LOW,  axis3) <= coordinates(i,j,k,axis3)) &
              .AND. (coordinates(i,j,k,axis3) <= gr_globalDomain(HIGH, axis3))) then
            regionData(i, j, k, var) = 0.0
          else
            regionData(i, j, k, var) = 1.1 * var
          end if
        end do
      end do
    end do

    deallocate(coordinates)
  end if

  ! Write data to points outside of domain and its boundary
  if      ((face == LOW) .AND. containsBoundary(face, axis)) then
    do     var = 1, regionSize(STRUCTSIZE)
      do     k = 1, regionSize(THIRD_DIR)
        do   j = 1, regionSize(SECOND_DIR)
          do i = 1, guard
            regionData(i, j, k, var) = 1.1 * var
          end do
        end do
      end do
    end do
  else if ((face == HIGH) .AND. containsBoundary(face, axis)) then
    do     var = 1, regionSize(STRUCTSIZE)
      do     k = 1, regionSize(THIRD_DIR)
        do   j = 1, regionSize(SECOND_DIR)
          do i = regionSize(BC_DIR)-NGUARD+1, regionSize(BC_DIR)
            regionData(i, j, k, var) = 1.1 * var
          end do
        end do
      end do
    end do
  end if

  applied = .TRUE.
end subroutine Grid_bcApplyToRegionSpecialized

