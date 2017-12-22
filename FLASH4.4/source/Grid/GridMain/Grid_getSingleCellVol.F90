!!****if* source/Grid/GridMain/Grid_getSingleCellVol
!!
!! NAME
!!  Grid_getSingleCellVol
!!
!! SYNOPSIS
!!
!!  Grid_getSingleCellVol(integer(IN) :: blockid,
!!                        integer(IN) :: beginCount,
!!                        integer(IN) :: point(MDIM), 
!!                        real(OUT)   :: cellVolume)
!!  
!! DESCRIPTION 
!!  
!!  Gets cell volumes for a single cell in a given block. 
!!
!!  
!! ARGUMENTS 
!!
!!  blockid - integer local blockid
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified,
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified,
!!               the first index, 1, is the leftmost guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!! 
!!  point(MDIM):
!!           specifies the point to return
!!   
!!           point(1) = i
!!           point(2) = j
!!           point(3) = k
!!
!!           If a problem is only 2d, point(3) is ignored.  For 1d problems
!!           point(2) and point(3) are ignored.
!!
!!
!!  cellVolume - real value containing the cell volume
!!
!!
!! NOTES
!! 
!!  Current implementations of this interface assume that all cells in a 
!!  dimension of a block have the same grid spacing. The grid spacings used
!!  are the ones returned by Grid_getDeltas.
!!
!! SEE ALSO
!!
!!  Grid_getDeltas
!!  Grid_getSingleCellCoords
!!
!!***

subroutine Grid_getSingleCellVol(blockID, beginCount, point, cellvolume)

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getSingleCellCoords, &
                             Grid_getBlkRefineLevel,                   &
                             Grid_getGeometry

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID, beginCount
  integer, intent(in) :: point(MDIM)
  real, intent(out) :: cellvolume

  integer :: geometry
  integer :: level
  real :: del(MDIM)
  real :: centerCoords(MDIM), leftCoords(MDIM), rightCoords(MDIM)
  
  del = 1.0
 
  call Grid_getGeometry(geometry)
  call Grid_getBlkRefineLevel(blockID, level)
  call Grid_getDeltas(level, del)

  select case (geometry)

  case (CARTESIAN)
     if(NDIM == 1) then
        cellvolume = del(IAXIS)
     else if(NDIM == 2) then
        cellvolume = del(IAXIS) * del(JAXIS)
     else
        cellvolume = del(IAXIS) * del(JAXIS) * del(KAXIS)
     end if

  case (POLAR)
     call Grid_getSingleCellCoords(point, blockID, CENTER, beginCount, centerCoords)

     if(NDIM == 1) then
        cellvolume = del(IAXIS) * 2.*PI * centerCoords(IAXIS)
     else if(NDIM == 2) then
        cellvolume = del(IAXIS) * del(JAXIS) * centerCoords(IAXIS)
     else
        cellvolume = del(IAXIS) * del(JAXIS) * centerCoords(IAXIS) * del(KAXIS)
     end if

  case (CYLINDRICAL)
     call Grid_getSingleCellCoords(point, blockID, CENTER, beginCount, centerCoords)

     if(NDIM == 1) then
        cellvolume = del(IAXIS) * 2.*PI * centerCoords(IAXIS)
     else if(NDIM == 2) then
        cellvolume = del(IAXIS) * 2.*PI * centerCoords(IAXIS) * del(JAXIS)
     else
        cellvolume = del(IAXIS) * del(JAXIS) * centerCoords(IAXIS) * del(KAXIS)
     end if

  case (SPHERICAL)
     call Grid_getSingleCellCoords(point, blockID, LEFT_EDGE, beginCount, leftCoords)
     call Grid_getSingleCellCoords(point, blockID, RIGHT_EDGE, beginCount, rightCoords)

     cellvolume = del(IAXIS) *  &
          ( leftCoords(IAXIS)*  leftCoords(IAXIS)  +  &
            leftCoords(IAXIS)* rightCoords(IAXIS)  +  &
           rightCoords(IAXIS)* rightCoords(IAXIS) )
     if(NDIM == 1) then
        cellvolume = cellvolume * 4.*PI/3.
     else if(NDIM == 2) then
        cellvolume = cellvolume * ( cos(leftCoords(JAXIS)) - cos(rightCoords(JAXIS)) ) * 2.*PI/3.
     else
        cellvolume = cellvolume * ( cos(leftCoords(JAXIS)) - cos(rightCoords(JAXIS)) ) *  &
             del(KAXIS) / 3.0
     end if

  end select

  return
end subroutine Grid_getSingleCellVol

subroutine Grid_getSingleCellVol_Itor(blockDesc, point, cellvolume, indexing)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getSingleCellCoords, &
                             Grid_getGeometry
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"

  type(block_metadata_t), intent(in) :: blockDesc
  integer, intent(in) :: point(MDIM)
  real, intent(out) :: cellvolume
  integer, intent(in),OPTIONAL :: indexing

  integer :: geometry
  integer :: blockID, level
  integer :: beginCount
  real    :: del(MDIM)
  real    :: centerCoords(MDIM), leftCoords(MDIM), rightCoords(MDIM)
 
  if (present(indexing)) then
     beginCount = indexing
  else
     beginCount = DEFAULTIDX
  end if

  blockID = blockDesc%id
  level = blockDesc%level
  if (level < 1 .OR. beginCount==EXTERIOR .OR. beginCount==INTERIOR) then
     call Grid_getSingleCellVol(blockID, beginCount, point, cellvolume)
     return
  end if

  if (beginCount .NE. GLOBALIDX1 .AND. beginCount .NE. DEFAULTIDX) then
     call Driver_abortFlash('Grid_getSingleCellVol_Itor: beginCount other than EXTERIOR not supported!')
  end if

  del = 1.0
  
  call Grid_getGeometry(geometry)
  call Grid_getDeltas(blockDesc%level, del)

  select case (geometry)

  case (CARTESIAN)
     if(NDIM == 1) then
        cellvolume = del(IAXIS)
     else if(NDIM == 2) then
        cellvolume = del(IAXIS) * del(JAXIS)
     else
        cellvolume = del(IAXIS) * del(JAXIS) * del(KAXIS)
     end if

  case (POLAR)
     if (blockID<1) call Driver_abortFlash("Grid_getSingleCellVol: got blockDesc without valid id")
     call Grid_getSingleCellCoords(point, level, CENTER, centerCoords)

     if(NDIM == 1) then
        cellvolume = del(IAXIS) * 2.*PI * centerCoords(IAXIS)
     else if(NDIM == 2) then
        cellvolume = del(IAXIS) * del(JAXIS) * centerCoords(IAXIS)
     else
        cellvolume = del(IAXIS) * del(JAXIS) * centerCoords(IAXIS) * del(KAXIS)
     end if

  case (CYLINDRICAL)
     call Grid_getSingleCellCoords(point, blockDesc, CENTER, beginCount, centerCoords)

     if(NDIM == 1) then
        cellvolume = del(IAXIS) * 2.*PI * centerCoords(IAXIS)
     else if(NDIM == 2) then
        cellvolume = del(IAXIS) * 2.*PI * centerCoords(IAXIS) * del(JAXIS)
     else
        cellvolume = del(IAXIS) * del(JAXIS) * centerCoords(IAXIS) * del(KAXIS)
     end if

  case (SPHERICAL)
     if (blockID<1) call Driver_abortFlash("Grid_getSingleCellVol: got blockDesc without valid id")
     call Grid_getSingleCellCoords(point, level, LEFT_EDGE, leftCoords)
     call Grid_getSingleCellCoords(point, level, RIGHT_EDGE, rightCoords)

     cellvolume = del(IAXIS) *  &
          ( leftCoords(IAXIS)*  leftCoords(IAXIS)  +  &
            leftCoords(IAXIS)* rightCoords(IAXIS)  +  &
           rightCoords(IAXIS)* rightCoords(IAXIS) )
     if(NDIM == 1) then
        cellvolume = cellvolume * 4.*PI/3.
     else if(NDIM == 2) then
        cellvolume = cellvolume * ( cos(leftCoords(JAXIS)) - cos(rightCoords(JAXIS)) ) * 2.*PI/3.
     else
        cellvolume = cellvolume * ( cos(leftCoords(JAXIS)) - cos(rightCoords(JAXIS)) ) *  &
             del(KAXIS) / 3.0
     end if

  end select

  return
end subroutine Grid_getSingleCellVol_Itor

