!!****if* source/Grid/GridMain/AMR/Amrex/gr_updateData
!!
!! NAME
!!
!!  gr_updateData
!!
!! SYNOPSIS
!!
!!  call gr_updateData()
!!
!! DESCRIPTION
!!
!!
!! This is a subroutine private to the Grid unit  which should be called
!! to updated some arrays that provide a rudimentary simulacrum of
!! same data  private to the PARAMESH implementation.
!! Sould be called  immediately after the mesh package finishes refining or
!! updating refinement.
!! It should be called only when there is re-gridding going on. For each
!! block on a processor, this routine gets the bounding box and refinement level
!! information, and some similar stuff.
!!
!! More specifically, a table of PARAMESHy metainformation arrays and potential Amrex-y
!! routines to get equivalent information:
!!
!!   lrefine     -  (Grid_getBlkRefineLevel unimplemented for desc)
!!   nodetype    -  (Grid_getBlkType unimplemented for desc, returns always LEAF for blockID)
!!   coord       -  Grid_getBlkCenterCoords (unimpl!), see Grid_getBlkBoundBox
!!   bsize       -  Grid_getBlkPhysicalSize (unimpl!), see Grid_getDeltas,Grid_getBlkBoundBox
!!   bnd_box     -  Grid_getBlkBoundBox[_desc]
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_updateData()
  use Grid_interface,        ONLY : Grid_getTileIterator, &
                                    Grid_releaseTileIterator
  use flash_iterator,        ONLY : flash_iterator_t
  use flash_tile,            ONLY : flash_tile_t
  use amrex_geometry_module, ONLY : amrex_problo
  use amrex_box_module,      ONLY : amrex_box
  use amrex_boxarray_module, ONLY : amrex_boxarray
  use amrex_amrcore_module,  ONLY : amrex_max_level, amrex_ref_ratio, amrex_get_boxarray

  use gr_specificData, ONLY : gr_ioLocalNumBlocks
  use gr_specificData, ONLY : gr_ioBlkNodeType, gr_ioBlkCoords, gr_ioBlkBsize , gr_ioBlkBoundBox, &
                              gr_ioBlkLrefine
  use Grid_data, ONLY: gr_meshMe, gr_boxContainingLeafNodes

  implicit none
  
  integer :: nBlocks, i
  real :: bnd_box_x,bnd_box_y,bnd_box_z

  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc
  type(amrex_box)        :: bx
  type(amrex_boxarray)   :: fba
  integer :: lb, lev
  integer :: rr
  logical :: blockHasChildren

  call Grid_getLocalNumBlks(nBlocks)
  gr_ioLocalNumBlocks = nBlocks

  if(allocated(gr_ioBlkLrefine) ) deallocate(gr_ioBlkLrefine)
  if(allocated(gr_ioBlkNodeType)) deallocate(gr_ioBlkNodeType)
  if(allocated(gr_ioBlkBoundBox)) deallocate(gr_ioBlkBoundBox)
  if(allocated(gr_ioBlkBsize)   ) deallocate(gr_ioBlkBsize)
  if(allocated(gr_ioBlkCoords)  ) deallocate(gr_ioBlkCoords)

  allocate(gr_ioBlkLrefine     (nBlocks))
  allocate(gr_ioBlkNodeType    (nBlocks))
  allocate(gr_ioBlkBoundBox    (LOW:HIGH,MDIM,nBlocks))
  allocate(gr_ioBlkBsize       (         MDIM,nBlocks))
  allocate(gr_ioBlkCoords      (         MDIM,nBlocks))

  call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)  ; lb = 1
  do while (itor%isValid())
     call itor%currentTile(tileDesc)
     gr_ioBlkNodeType(lb) = LEAF ! DEV: NOT ALWAYS TRUE!!!
     gr_ioBlkLrefine(lb)  = tileDesc % level
     lev                  = tileDesc % level - 1
     if (lev .GE. amrex_max_level) then
        blockHasChildren = .FALSE.
     else
        fba = amrex_get_boxarray(lev+1)
        rr = amrex_ref_ratio(lev)
        bx = amrex_box(tileDesc % limits(LOW,1:NDIM )-1, &
                       tileDesc % limits(HIGH,1:NDIM)-1)
        call bx%refine(rr)   !Note: this modifies bx, do not use naively after this!
        blockHasChildren = fba%intersects(bx)
     end if
     if (blockHasChildren) then
        gr_ioBlkNodeType(lb) = PARENT_BLK ! DEV: Do we need to distinguish ANCESTOR?
     else
        gr_ioBlkNodeType(lb) = LEAF
     end if

     call tileDesc%boundBox(gr_ioBlkBoundBox(:,:,lb))
     gr_ioBlkBsize(:,lb) =  gr_ioBlkBoundBox(HIGH,:,lb) - gr_ioBlkBoundBox(LOW,:,lb)
     gr_ioBlkCoords(:,lb)= (gr_ioBlkBoundBox(HIGH,:,lb)+gr_ioBlkBoundBox(LOW,:,lb))*0.5

     call itor%next()              ; lb = lb+1
  enddo
  call Grid_releaseTileIterator(itor)

  nBlocks = lb-1
  if (gr_ioLocalNumBlocks .NE. nBlocks) then
     print*,'gr_ioLocalNumBlocks,nBlocks=',gr_ioLocalNumBlocks,nBlocks,' @',gr_meshMe
     call Driver_abortFlash("gr_ioLocalNumBlocks .NE. nBlocks")
    end if



  gr_boxContainingLeafNodes(LOW ,IAXIS) =  HUGE(bnd_box_x)
  gr_boxContainingLeafNodes(HIGH,IAXIS) = -HUGE(bnd_box_x)
#if NDIM > 1
  gr_boxContainingLeafNodes(LOW ,JAXIS) =  HUGE(bnd_box_x)
  gr_boxContainingLeafNodes(HIGH,JAXIS) = -HUGE(bnd_box_x)
#else
  gr_boxContainingLeafNodes(LOW ,JAXIS) =  amrex_problo(JAXIS)
  gr_boxContainingLeafNodes(HIGH,JAXIS) =  amrex_problo(JAXIS)
#endif
#if NDIM > 2
  gr_boxContainingLeafNodes(LOW ,KAXIS) =  HUGE(bnd_box_x)
  gr_boxContainingLeafNodes(HIGH,KAXIS) = -HUGE(bnd_box_x)
#else
  gr_boxContainingLeafNodes(LOW ,KAXIS) =  amrex_problo(KAXIS)
  gr_boxContainingLeafNodes(HIGH,KAXIS) =  amrex_problo(KAXIS)
#endif

  do i = 1,nBlocks
     bnd_box_x = gr_ioBlkBoundBox(1,1,i)
     bnd_box_y = gr_ioBlkBoundBox(1,2,i)
     bnd_box_z = gr_ioBlkBoundBox(1,3,i)
     if (gr_ioBlkNodeType(i) == 1) then ! for LEAF blocks
        gr_boxContainingLeafNodes(LOW ,IAXIS) = min(bnd_box_x,              gr_boxContainingLeafNodes(LOW,IAXIS))
        gr_boxContainingLeafNodes(HIGH,IAXIS) = max(gr_ioBlkBoundBox(2,1,i),gr_boxContainingLeafNodes(HIGH,IAXIS))
#if NDIM > 1
        gr_boxContainingLeafNodes(LOW ,JAXIS) = min(bnd_box_y,              gr_boxContainingLeafNodes(LOW,JAXIS))
        gr_boxContainingLeafNodes(HIGH,JAXIS) = max(gr_ioBlkBoundBox(2,2,i),gr_boxContainingLeafNodes(HIGH,JAXIS))
#endif
#if NDIM > 2
        gr_boxContainingLeafNodes(LOW ,KAXIS) = min(bnd_box_z,              gr_boxContainingLeafNodes(LOW,KAXIS))
        gr_boxContainingLeafNodes(HIGH,KAXIS) = max(gr_ioBlkBoundBox(2,3,i),gr_boxContainingLeafNodes(HIGH,KAXIS))
#endif
     end if

  end do

end subroutine gr_updateData


