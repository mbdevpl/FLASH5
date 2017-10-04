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
    !!  block_iterator_t itor = block_iterator_t(integer(IN)         :: nodetype,
    !!                                           level(IN), optional :: level)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks
    !!  within the current paramesh octree structure.  The iterator is already
    !!  set to the first matching block.
    !!
    !! ARGUMENTS
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS)
    !!  level    - if nodetype is LEAF, PARENT, ANCESTOR, or REFINEMENT, then 
    !!             iterate only over blocks located at this level of 
    !!             octree structure
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    function init_iterator(nodetype, level) result(this)
        integer, intent(IN)           :: nodetype
        integer, intent(IN), optional :: level
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

