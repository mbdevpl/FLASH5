!!****f* source/Grid/GridMain/AMR/Amrex/Grid_getBlkBoundBox
!!
!! NAME
!!  Grid_getBlkBoundBox
!!
!! SYNOPSIS
!!
!! 
!!  Grid_getBlkBoundBox(integer(IN)  :: blockDesc,
!!                         real(OUT) :: boundBox(2, MDIM))
!!  
!! DESCRIPTION 
!!
!!  Gets the physical domain bounding box of the block identified 
!!  by blockId.  For each dimension the left (lower or forward) 
!!  physical coordinate of the block edge and the right (upper or back) 
!!  physical coordinate of the block edge is returned.  See arguments
!!  below for more detail.
!!
!! ARGUMENTS
!!
!!  blockId -local block number
!!
!!  boundBox - returned array holding the boundBox coordinates in
!!             each dimension
!!
!!            for readability, in constants.h we define IAXIS = 1, JAXIS = 2, KAXIS = 3
!!
!!            boundBox(1,IAXIS) = left edge coordinate of block in x direction
!!            boundBox(2,IAXIS) = right edge coordinate of block in x direction
!!            boundBox(1,JAXIS) = top edge coordinate of block in y direction
!!            boundBox(2,JAXIS) = bottom edge coordinate of block in y direction
!!            boundBox(1,KAXIS) = front edge coordinate of block in z direction
!!            boundBox(2,KAXIS) = back edge coordinate of block in z direction
!!
!! EXAMPLE
!!  
!!   In 2 dimensions, if physical coordinates are ...
!!    
!!     ________________(0.5 1.0)
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |_______________ |
!!  (-0.5, 0.0)
!!
!!
!!
!!     boundBox(1, IAXIS) = -0.5
!!     boundBox(2, IAXIS) = 0.5
!!     boundBox(1, JAXIS) = 0.0
!!     boundBox(2, JAXIS) = 1.0
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkBoundBox_desc(blockDesc, boundBox)
  use amrex_amrcore_module,  ONLY : amrex_geom
  use amrex_geometry_module, ONLY : amrex_problo

  use block_metadata,        ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(IN)  :: blockDesc
  real,                   intent(OUT) :: boundBox(LOW:HIGH, MDIM)

  ! DEV: FIXME How to manage matching amrex_real to FLASH real
  boundBox = 1.0
  associate(x0   => amrex_problo, &
            dx   => amrex_geom(blockDesc%level - 1)%dx, &
            lo   => blockDesc%limits(LOW,  :), &
            hi   => blockDesc%limits(HIGH, :))
    ! lo is 1-based cell-index of lower-left cell in block 
    ! hi is 1-based cell-index of upper-right cell in block
    boundBox(LOW,  1:NDIM) = x0(1:NDIM) + (lo(1:NDIM) - 1)*dx(1:NDIM)
    boundBox(HIGH, 1:NDIM) = x0(1:NDIM) + (hi(1:NDIM)    )*dx(1:NDIM)
  end associate
end subroutine Grid_getBlkBoundBox_desc

