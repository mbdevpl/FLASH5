!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFaceVar/Grid_bcApplyToRegionSpecialized
!!
!! NAME
!!  Grid_bcApplyToRegionSpecialized
!!
!! SYNOPSIS
!!  call Grid_bcApplyToRegionSpecialized(integer(IN)          :: bcType, 
!!                                       integer(IN)          :: gridDataStruct,
!!                                       integer(IN)          :: guard,
!!                                       integer(IN)          :: axis,
!!                                       integer(IN)          :: face,
!!                                       real(INOUT)          :: regionData(regionSize(BC_DIR),
!!                                                                          regionSize(SECOND_DIR), 
!!                                                                          regionSize(THIRD_DIR),
!!                                                                          regionSize(STRUCTSIZE))
!!                                       integer(IN)          :: regionSize(REGION_DIM),
!!                                       logical(IN)          :: mask(STRUCTSIZE),
!!                                       logical(OUT)         :: applied,
!!                                       block_metadata_t(IN) :: blockDesc,
!!                                       integer(IN)          :: secondDir,
!!                                       integer(IN)          :: thirdDir,
!!                                       integer(IN)          :: endPoints(LOW:HIGH, MDIM),
!!                             optional, integer(IN)          :: idest)
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

subroutine Grid_bcApplyToRegionSpecialized(bcType, gridDataStruct, &
                                           guard, axis, face, &
                                           regionData, regionSize, &
                                           mask, applied, &
                                           blockDesc, &
                                           secondDir, thirdDir, &
                                           endPoints, idest)
  use Grid_data,        ONLY : gr_globalDomain 
  use Grid_interface,   ONLY : Grid_getBlkBC, &
                               Grid_getSingleCellCoords, &
                               Grid_getBlkBoundBox, &
                               Grid_getDeltas
  use block_metadata,   ONLY : block_metadata_t

  implicit none

  integer,                intent(IN)           :: bcType
  integer,                intent(IN)           :: axis
  integer,                intent(IN)           :: face
  integer,                intent(IN)           :: guard
  integer,                intent(IN)           :: gridDataStruct
  integer,                intent(IN)           :: regionSize(REGION_DIM)
  real,                   intent(INOUT)        :: regionData(regionSize(BC_DIR),     &
                                                             regionSize(SECOND_DIR), &
                                                             regionSize(THIRD_DIR),  &
                                                             regionSize(STRUCTSIZE))
  logical,                intent(IN)           :: mask(regionSize(STRUCTSIZE))
  logical,                intent(OUT)          :: applied
  type(block_metadata_t), intent(IN)           :: blockDesc
  integer,                intent(IN)           :: secondDir,thirdDir
  integer,                intent(IN)           :: endPoints(LOW:HIGH, MDIM)
  integer,                intent(IN), OPTIONAL :: idest

  integer :: i, j, k, var
  integer :: axis2, axis3
  logical :: isOnBoundary
  logical :: containsBoundary(LOW:HIGH, 1:MDIM)

  real    :: x0(1:MDIM)
  real    :: blkBounds(LOW:HIGH, 1:MDIM)
  real    :: deltas(1:MDIM)
  real    :: coord(1:MDIM)

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

  write(*,*) "Box Low:       ", blockDesc%limits(LOW, :)
  write(*,*) "Box High:      ", blockDesc%limits(HIGH, :)
#endif

  call Grid_getBlkBoundBox(blockDesc, blkBounds)
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
    call Grid_getSingleCellCoords([1,1,1], blockDesc, CENTER, INTERIOR, x0)
    call Grid_getDeltas(blockDesc%level, deltas)

    ! Set all face-centered values on domain boundary to zero.  Those
    ! coplanar to, but not on the boundary are considered guardcells
    if (face == LOW) then
      i = guard+1
    else
      i = regionSize(BC_DIR) - NGUARD
    end if

    do     var = 1, regionSize(STRUCTSIZE)
      do     k = 1, regionSize(THIRD_DIR)
        do   j = 1, regionSize(SECOND_DIR)
          ! Determine location of current face in global index space
          if      (axis == IAXIS) then
            coord(IAXIS) = x0(IAXIS) + (i - 1.5)*deltas(IAXIS)
            coord(JAXIS) = x0(JAXIS) + (j - 1.0)*deltas(JAXIS)
            coord(KAXIS) = x0(KAXIS) + (k - 1.0)*deltas(KAXIS)
            axis2 = JAXIS
            axis3 = KAXIS
          else if (axis == JAXIS) then
            coord(IAXIS) = x0(IAXIS) + (j - 1.0)*deltas(IAXIS)
            coord(JAXIS) = x0(JAXIS) + (i - 1.5)*deltas(JAXIS)
            coord(KAXIS) = x0(KAXIS) + (k - 1.0)*deltas(KAXIS)
            axis2 = IAXIS
            axis3 = KAXIS
          else
            coord(IAXIS) = x0(IAXIS) + (j - 1.0)*deltas(IAXIS)
            coord(JAXIS) = x0(JAXIS) + (k - 1.0)*deltas(JAXIS)
            coord(KAXIS) = x0(KAXIS) + (i - 1.5)*deltas(KAXIS)
            axis2 = IAXIS
            axis3 = JAXIS
          end if

          ! By construction, the face is coplanar to boundary.  Check if
          ! it is actually on boundary.
          if (      (gr_globalDomain(LOW,  axis2) <= coord(axis2)) &
              .AND. (coord(axis2) <= gr_globalDomain(HIGH, axis2)) &
              .AND. (gr_globalDomain(LOW,  axis3) <= coord(axis3)) &
              .AND. (coord(axis3) <= gr_globalDomain(HIGH, axis3))) then
            regionData(i, j, k, var) = 0.0
          else
            regionData(i, j, k, var) = 1.1 * var
          end if
        end do
      end do
    end do
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

