!!****ih* source/Grid/GridMain/AMR/Amrex/block_1lev_iterator
!!
!! This module is a facade pattern that maps the AMReX Fortran iterator onto 
!! the interface required presently by FLASH.
!!
!! Ideally, we will be able to use the AMReX iterator directly in the code 
!! and client code will gain access to it through implementation-specific 
!! code like Grid_getBlkIterator.
!!
!! This is a variant that uses type(amrex_mfiter) as underlying iterator.
!! It iterates only over block at the same, given refinement level.
!!****

#include "constants.h"

module block_1lev_iterator

    use amrex_multifab_module, ONLY : amrex_mfiter, &
                                    amrex_mfiter_build, &
                                    amrex_mfiter_destroy

    implicit none

    private

    !!****ic* block_1lev_iterator/block_1lev_iterator_t
    !!
    !! NAME
    !!  block_1lev_iterator_t
    !!
    !!****
    type, public :: block_1lev_iterator_t
        type(amrex_mfiter),POINTER :: mfi   => NULL()
        integer                 :: level    = INVALID_LEVEL
        logical                 :: isValid = .FALSE.
    contains
        procedure, public :: is_valid
        procedure, public :: first
        procedure, public :: next
        procedure, public :: tilebox
        procedure, public :: blkMetaData
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
        final             :: destroy_iterator
#endif
    end type block_1lev_iterator_t

    interface block_1lev_iterator_t
        procedure :: init_iterator_mf
        procedure :: init_iterator
        procedure :: init_iterator_mfa
    end interface block_1lev_iterator_t

contains

    !!****im* block_1lev_iterator_t/block_1lev_iterator_t
    !!
    !! NAME
    !!  block_1lev_iterator_t
    !!
    !! SYNOPOSIS
    !!  block_1lev_iterator_t itor = block_1lev_iterator_t(integer(IN)         :: nodetype,
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
  function init_iterator_mf(nodetype, mf, level, tiling) result(this)
    use amrex_multifab_module, ONLY : amrex_multifab

        type(block_1lev_iterator_t)        :: this
        integer, intent(IN)           :: nodetype
        type(amrex_multifab),intent(IN) :: mf
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling

        if (present(level)) then
            this%level = level
        end if
 
        allocate(this%mfi)

        ! DEVNOTE: the AMReX iterator is not built based on nodetype.
        ! It appears that we get leaves every time.

        ! Initial iterator is not primed.  Advance to first compatible block.
        call amrex_mfiter_build(this%mfi,mf,tiling=tiling)
        call this%next()
  end function init_iterator_mf

  function init_iterator_mfa(nodetype, mfArray, level, tiling) result(this)
    use amrex_multifab_module, ONLY : amrex_multifab

        type(block_1lev_iterator_t)        :: this
        integer, intent(IN)           :: nodetype
        type(amrex_multifab),intent(IN) :: mfArray(*)
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling

        if (present(level)) then
            this%level = level
        end if
 
        allocate(this%mfi)

        ! DEVNOTE: the AMReX iterator is not built based on nodetype.
        ! It appears that we get leaves every time.

        ! Initial iterator is not primed.  Advance to first compatible block.
        call amrex_mfiter_build(this%mfi,mfArray(level),tiling=tiling)
        call this%next()
     end function init_iterator_mfa

    function init_iterator(nodetype, level, tiling) result(this)
      use amrex_multifab_module, ONLY : amrex_multifab
      use gr_amrextData

        type(block_1lev_iterator_t)        :: this
        integer, intent(IN)           :: nodetype
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling

        type(amrex_multifab),POINTER :: mfArray(:)
        type(amrex_multifab),POINTER :: mf

        if (present(level)) then
            this%level = level
        end if
 
        allocate(this%mfi)

!!$        call gr_amrextDataInit(10)
!!$        call gr_amrextBuildMultiFabsFromF4Grid(gr_amrextUnkMFs,lrefine_max,LEAF)
        if (present(level)) then
           mfArray => gr_amrextUnkMFs
           call amrex_mfiter_build(this%mfi,mfArray(level),tiling=tiling)
        else
           mf => gr_amrextUnkMFs(1)
           call amrex_mfiter_build(this%mfi,mf,tiling=tiling)
        end if
        call this%next()
    end function init_iterator

#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
    !!****im* block_1lev_iterator_t/destroy_iterator
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
        type(block_1lev_iterator_t), intent(INOUT) :: this

        call amrex_mfiter_destroy(this%mfi)
        deallocate(this%mfi)
        nullify(this%mfi)
    end subroutine destroy_iterator
#endif

    !!****m* block_1lev_iterator_t/first
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
        class(block_1lev_iterator_t), intent(INOUT) :: this

        ! reset to before first valid block
        call this%mfi%clear()
        this%isValid = .FALSE.  !DEV: ??
        ! Initial iterator is not primed.  Advance to first compatible block.
        call this%next()
    end subroutine first
 
    !!****m* block_1lev_iterator_t/is_valid
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
    function is_valid(this) result(ans)
        class(block_1lev_iterator_t), intent(IN) :: this
        logical :: ans

        ans = this%isValid
    end function is_valid

    !!****m* block_1lev_iterator_t/next
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
        class(block_1lev_iterator_t), intent(INOUT) :: this

        this%isValid = this%mfi%next()

    end subroutine next

    function tilebox (this) result (bx)
      use amrex_box_module, ONLY : amrex_box
      class(block_1lev_iterator_t), intent(in) :: this
      type(amrex_box) :: bx
      integer :: inodal(3)
      inodal = 0
      bx = this%mfi%tilebox()
    end function tilebox

    !!****m* block_1lev_iterator_t/blkMetaData
    !!
    !! NAME
    !!  blkMetaData 
    !!
    !! SYNPOSIS
    !!  call itor%blkMetaData(block_metadata_t(OUT) : blockDesc)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block currently set in the
    !!  iterator.
    !!
    !!****
    subroutine blkMetaData(this, blockDesc)
        use amrex_box_module, ONLY : amrex_box
        use block_metadata, ONLY : block_metadata_t

        class(block_1lev_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: blockDesc

        type(amrex_box) :: box
       
        box = this%mfi%tilebox()

        ! TODO: Determine if box contains GC or not and finalize limits/limitsGC
!!$        blockDesc%grid_index        = this%oti%grid_index()
        blockDesc%level             = this%level
        blockDesc%limits(LOW, :)    = box%lo
        blockDesc%limits(HIGH, :)   = box%hi
        blockDesc%limitsGC(LOW, :)  = box%lo
        blockDesc%limitsGC(HIGH, :) = box%hi
    end subroutine blkMetaData
 
end module block_1lev_iterator

