!!****ih* source/Grid/gr_iterator
!!
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"

module gr_iterator

    implicit none

    private

    public :: build_iterator, destroy_iterator

    !!****ic* gr_iterator/leaf_iterator_t
    !!
    !! NAME
    !!  gr_iterator_t
    !!
    !!****
    type, public :: gr_iterator_t
    contains
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: blkMetaData
    end type gr_iterator_t

contains

    !!****im* gr_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(gr_iterator_t(OUT) :: itor,
    !!                 integer(IN), optional :: level,
    !!                 logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of leaf blocks or
    !!  tiles within the current octree structure.  The iterator is already
    !!  set to the first matching leaf block/tile.
    !!
    !! ARGUMENTS
    !!  itor     - the constructed iterator
    !!  level    - iterate only over leaf blocks/tiles located at this level of
    !!             refinement.
    !!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
    !!             walk across all associated blocks on a tile-by-tile basis *if*
    !!             the implementation supports this feature.  If a value is not
    !!             given, is FALSE, or the implementation does not support tiling,
    !!             the iterator will iterate on a block-by-block basis.
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine build_iterator(itor, level, tiling)
        use Driver_interface, ONLY : Driver_abortFlash

        type(gr_iterator_t), intent(OUT)          :: itor
        integer,               intent(IN), optional :: level
        logical,               intent(IN), optional :: tiling
        
        call Driver_abortFlash("[build_iterator] You are working with a useless gr_iterator_t stub")
    end subroutine build_iterator

    !!****im* gr_iterator_t/destroy_iterator
    !!
    !! NAME
    !!  destroy_iterator
    !!
    !! SYNOPOSIS
    !!  destroy_iterator(gr_iterator_t(INOUT) :: itor)
    !!
    !! DESCRIPTION
    !!  Destroy the given iterator.
    !!
    !! ARGUMENTS
    !!  itor     - the iterator to destroy
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(itor)
        type(gr_iterator_t), intent(INOUT) :: itor

    end subroutine destroy_iterator

    !!****m* gr_iterator_t/is_valid
    !!
    !! NAME
    !!  is_valid
    !!
    !! SYNPOSIS
    !!  logical valid = itor%is_valid()
    !!
    !! DESCRIPTION
    !!  Determine if the iterator is currently set to a valid block.
    !!
    !! RETURN VALUE 
    !!  True if iterator is currently set to a valid block
    !!
    !!****
    logical function is_valid(this)
        class(gr_iterator_t), intent(IN) :: this

        write(*,*) "[is_valid] You are working with a useless gr_iterator_t stub"
        stop
    end function is_valid

    !!****m* gr_iterator_t/next
    !!
    !! NAME
    !!  next
    !!
    !! SYNPOSIS
    !!  call itor%next()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next block managed by process and that meets
    !!  the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next(this)
        class(gr_iterator_t), intent(INOUT) :: this

        write(*,*) "[next] You are working with a useless gr_iterator_t stub"
        stop
    end subroutine next

    !!****m* gr_iterator_t/blkMetaData
    !!
    !! NAME
    !!  blkMetaData 
    !!
    !! SYNPOSIS
    !!  call itor%blkMetaData(block_metadata_t(OUT) : block)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block currently set in the
    !!  iterator.
    !!
    !!****
    subroutine blkMetaData(this, mData)
        use block_metadata, ONLY : block_metadata_t

        class(gr_iterator_t), intent(IN)  :: this
        type(block_metadata_t), intent(OUT) :: mData

        write(*,*) "[blkMetaData] You are working with a useless gr_iterator_t stub"
        stop
    end subroutine blkMetaData
    
  end module gr_iterator

