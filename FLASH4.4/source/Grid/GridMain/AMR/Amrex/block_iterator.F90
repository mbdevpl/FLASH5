!!****ih* source/Grid/GridMain/AMR/Amrex/block_iterator
!!
!! This module is a facade pattern that maps the AMReX Fortran iterator onto 
!! the interface required presently by FLASH.
!!
!! Ideally, we will be able to use the AMReX iterator directly in the code 
!! and client code will gain access to it through implementation-specific 
!! code like Grid_getBlkIterator.
!!
!!****

module block_iterator

    use amrex_octree_module, ONLY : amrex_octree_iter, &
                                    amrex_octree_iter_build, &
                                    amrex_octree_iter_destroy

    implicit none

    private

    !!****ic* block_iterator/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !!****
    type, public :: block_iterator_t
        type(amrex_octree_iter) :: oti
        integer                 :: level    = INVALID_LEVEL
        logical                 :: is_valid = .FALSE.
    contains
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: blkMetaData
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
        final             :: destroy_iterator
#endif
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

        if (present(level)) then
            this%level = level
        end if
 
        ! DEVNOTE: the AMReX iterator is not built based on nodetype.
        ! It appears that we get leaves every time.

        ! Initial iterator is not primed.  Advance to first compatible block.
        call amrex_octree_iter_build(this%oti)
        call this%next()
    end function init_iterator

#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
    !!****im* block_iterator_t/destroy_iterator
    !!
    !! NAME
    !!  destroy_iterator
    !!
    !! SYNPOSIS
    !!  Called automatically
    !!
    !! DESCRIPTION
    !!  Clean-up block interator object at destruction
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(this)
        type(block_iterator_t), intent(INOUT) :: this

        call amrex_octree_iter_destroy(this%oti)
    end subroutine destroy_iterator
#endif

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

        return this%is_valid
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

        this%is_valid = this%oti%next()

        if (this%level /= INVALID_LEVEL) then
            ! Search for leaves on given level
            do while (this%is_valid)
                if (this%oti%level() == this%level) then
                    exit
                end if

                this%is_valid = this%oti%next()
            end do
        end if
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
    subroutine blkMetaData(this, block)
        use block_metadata, ONLY : block_metadata_t

        class(block_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: block

        type(amrex_box) :: box
       
        box = this%oti%box()

        ! TODO: Determine if box contains GC or not and finalize limits/limitsGC
        block.grid_index        = this%oti%grid_index()
        block.level             = this%oti%level()
        block.limits(LOW, :)    = box%lo
        block.limits(HIGH, :)   = box%hi
        block.limitsGC(LOW, :)  = box%lo
        block.limitsGC(HIGH, :) = box%hi
    end subroutine blkMetaData
 
end module block_iterator

