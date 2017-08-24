subroutine gr_fillMetaData(blockID, blockDesc)
  use Grid_interface
  use block_metadata, ONLY : block_metadata_t
  implicit none
  integer,intent(IN) :: blockID
  type(block_metadata_t),intent(OUT) :: blockDesc

  blockDesc%id = blockID
  blockDesc%grid_index = -1 !DEV: ??
  call Grid_getBlkRefineLevel(blockID, blockDesc%level)
  call Grid_getBlkIndexLimits(blockID, blockDesc%limits, blockDesc%limitsGC)
  call Grid_getBlkCornerID   (blockID, blockDesc%cid,    blockDesc%stride)
end subroutine gr_fillMetaData
