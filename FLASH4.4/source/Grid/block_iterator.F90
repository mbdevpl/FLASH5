!!****ih* source/Grid/block_iterator
!!
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"

module block_iterator

    implicit none

    private

    public :: destroy_iterator

    !!****ic* block_iterator/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !!****
    type, public :: block_iterator_t
    contains
        procedure, public :: first
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: blkMetaData
    end type block_iterator_t

    interface block_iterator_t
        procedure :: init_iterator
    end interface block_iterator_t

contains

    !!****im* block_iterator_t/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !! SYNOPOSIS
    !!  block_iterator_t itor = block_iterator_t(integer(IN)           :: nodetype,
    !!                                           integer(IN), optional :: level,
    !!                                           logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the current octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!
    !! ARGUMENTS
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS)
    !!  level    - if nodetype is LEAF, PARENT, ANCESTOR, or REFINEMENT, then 
    !!             iterate only over blocks/tiles located at this level of
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
    function init_iterator(nodetype, level, tiling) result(this)
        integer, intent(IN)           :: nodetype
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling
        type(block_iterator_t)        :: this

        write(*,*) "You are working with a useless block_iterator_t stub"
        stop
    end function init_iterator

    !!****m* block_iterator_t/first
    !!
    !! NAME
    !!  first
    !!
    !! SYNPOSIS
    !!  call itor%first() 
    !!
    !! DESCRIPTION
    !!  Reset iterator to the initial block managed by process
    !!
    !!****
    subroutine first(this)
        class(block_iterator_t), intent(INOUT) :: this

        write(*,*) "You are working with a useless block_iterator_t stub"
        stop
    end subroutine first
 
    IMPURE_ELEMENTAL subroutine destroy_iterator(this)
        type(block_iterator_t), intent(INOUT) :: this

    end subroutine destroy_iterator

    !!****m* block_iterator_t/is_valid
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
        class(block_iterator_t), intent(IN) :: this

        write(*,*) "You are working with a useless block_iterator_t stub"
        stop
    end function is_valid

    !!****m* block_iterator_t/next
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
        class(block_iterator_t), intent(INOUT) :: this

        write(*,*) "You are working with a useless block_iterator_t stub"
        stop
    end subroutine next

    !!****m* block_iterator_t/blkMetaData
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
        class(block_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: mData

        write(*,*) "You are working with a useless block_iterator_t stub"
        stop
    end subroutine blkMetaData
    
end module block_iterator

