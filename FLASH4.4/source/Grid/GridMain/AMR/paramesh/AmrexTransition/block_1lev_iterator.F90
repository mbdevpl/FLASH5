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

#include "FortranLangFeatures.fh"
#include "constants.h"

module block_1lev_iterator

    use amrex_multifab_module, ONLY : amrex_multifab
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
        real,POINTER            :: fp(:,:,:,:)
        type(amrex_multifab),POINTER :: mf  => NULL()
    contains
        procedure, public :: is_valid
        procedure, public :: first
        procedure, public :: next
        procedure, public :: grid_index
        procedure, public :: tilebox
        procedure, public :: fabbox
        procedure, public :: dataPtr
        procedure, public :: blkMetaData
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
        final             :: destroy_iterator
#else
        procedure, public :: destroy_iterator
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
    !!             octree structure. !DEVNOTE: nodetype not implemented!
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
  function init_iterator_mf(nodetype, mf, level, tiling) result(this)
    use amrex_multifab_module, ONLY : amrex_multifab

        type(block_1lev_iterator_t)        :: this
        integer, intent(IN)           :: nodetype
        type(amrex_multifab),intent(IN),TARGET :: mf
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
        this%mf => mf
!!$        print*,'block_1lev_iterator: init_iterator_mf  on this=',this%isValid,this%level,associated(this%mfi)
        this%isValid = .TRUE.
        call this%next()
  end function init_iterator_mf

  function init_iterator_mfa(nodetype, mfArray, level, tiling) result(this)
    use amrex_multifab_module, ONLY : amrex_multifab

        type(block_1lev_iterator_t)        :: this
        integer, intent(IN)           :: nodetype
        type(amrex_multifab),intent(IN),TARGET :: mfArray(0:*)
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling

        if (present(level)) then
            this%level = level
        end if
 
        allocate(this%mfi)

        ! DEVNOTE: the AMReX iterator is not built based on nodetype.
        ! It appears that we get leaves every time.

        ! Initial iterator is not primed.  Advance to first compatible block.
        call amrex_mfiter_build(this%mfi,mfArray(level-1),tiling=tiling)
        this%mf => mfArray(level-1)
!!$        print*,'block_1lev_iterator: init_iterator_mfa on this=',this%isValid,this%level,associated(this%mfi)
        this%isValid = .TRUE.
        call this%next()
     end function init_iterator_mfa

    function init_iterator(nodetype, level, tiling) result(this)
      use amrex_multifab_module, ONLY : amrex_multifab
      use gr_physicalMultifabs,  ONLY : Unk

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

        if (present(level)) then
           mfArray => Unk
           mf => Unk(level-1)
!!$           print*,'amrex_multifab_nghost(mf)=',mf%nghost()
!!$           print*,'ABOUT TO call amrex_mfiter_build,size(mfArray)=',size(mfArray)
!!$           call amrex_mfiter_build(this%mfi,mfArray(level),tiling=tiling)
           call amrex_mfiter_build(this%mfi,mf,tiling=tiling)
!!$           print*,'amrex_multifab_nghost(mf)=',mf%nghost()
           this%mf => mfArray(level-1)
        else
           mf => Unk(0)
           call amrex_mfiter_build(this%mfi,mf,tiling=tiling)
           this%mf => mf
        end if
!!$        print*,'block_1lev_iterator: init_iterator     on this=',this%isValid,this%level,associated(this%mfi)
        this%isValid = .TRUE.
        call this%next()
    end function init_iterator

!#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
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
      class(block_1lev_iterator_t), intent(INOUT) :: this

      if (this%isValid) then
         call amrex_mfiter_destroy(this%mfi)
         deallocate(this%mfi)
         nullify(this%mfi)
         nullify(this%mf)
         this%isValid = .FALSE.
      end if
    end subroutine destroy_iterator
!#endif

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

        print*,'block_1lev_%first: IGNORING ATTEMPT  on this=',this%isValid,this%level,associated(this%mfi)

!!$        ! reset to before first valid block
!!$        print*,'block_1lev_%first: about to do clear on this=',this%isValid,this%level,associated(this%mfi)
!!$        call this%mfi%clear()
!!$        this%isValid = .TRUE.
!!$        ! Initial iterator is not primed.  Advance to first compatible block.
!!$        print*,'block_1lev_%first: about to do next on this=',this%isValid,this%level,associated(this%mfi)
!!$        call this%next()
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

        logical :: v
        if (this%isValid) then
!!$           print*,'block_1lev_iterator: about to do next on this=',this%isValid,this%level,associated(this%mfi)
           v = this%mfi%next()
!!$           print*,'block_1lev_iterator:       done  next on this=',       v    ,this%level,associated(this%mfi)

           if (.NOT. v) then
              call this%destroy_iterator()
           end if

        else
           print*,'block_1lev_iterator: no next, inValid! on this=',this%isValid,this%level,associated(this%mfi)
           call Driver_abortFlash("block_1lev_iterator: attempting next() on invalid!")
        end if

    end subroutine next

    function grid_index(this) result(idx)
      class(block_1lev_iterator_t), intent(IN) :: this
      integer :: idx
      idx = this%mfi%grid_index()
    end function grid_index

    function tilebox (this) result (bx)
      use amrex_box_module, ONLY : amrex_box
      class(block_1lev_iterator_t), intent(in) :: this
      type(amrex_box) :: bx
      integer :: inodal(3)
      inodal = 0
      bx = this%mfi%tilebox()
    end function tilebox

    function fabbox (this) result (bx)
      use amrex_box_module, ONLY : amrex_box
      class(block_1lev_iterator_t), intent(in) :: this
      type(amrex_box) :: bx
      integer :: inodal(3)
      inodal = 0
      bx = this%mfi%fabbox()
    end function fabbox

    function dataPtr (this) result (dp)
      use amrex_multifab_module, ONLY : amrex_multifab
      real, contiguous, pointer, dimension(:,:,:,:) :: dp
      class(block_1lev_iterator_t), intent(in) :: this
      dp => this%mf%dataPtr(this%mfi)
    end function dataPtr

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
        use tree,           ONLY : lrefine_max

        class(block_1lev_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: blockDesc

        type(amrex_box) :: box, fabbox
       
        box = this%mfi%tilebox()
        fabbox = this%mfi%fabbox()

        ! TODO: Determine if box contains GC or not and finalize limits/limitsGC
        blockDesc%grid_index        = this%mfi%grid_index()
        blockDesc%level             = this%level
        blockDesc%limits(LOW, :)    = box%lo
        blockDesc%limits(HIGH, :)   = box%hi
        blockDesc%limitsGC(LOW, :)  = fabbox%lo
        blockDesc%limitsGC(HIGH, :) = fabbox%hi

        blockDesc%localLimits(LOW, :)   = blockDesc%limits(LOW, :)   - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimits(HIGH, :)  = blockDesc%limits(HIGH, :)  - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimitsGC(LOW, :) = 1
        blockDesc%localLimitsGC(HIGH, :)= blockDesc%limitsGC(HIGH, :)- blockDesc%limitsGC(LOW, :) + 1

        blockDesc%fp => this%dataPtr()

        blockDesc%id = -888
        blockDesc%stride = 2**(lrefine_max - blockDesc%level)
        blockDesc%cid = (box%lo - 1) * blockDesc%stride + 1

    end subroutine blkMetaData
 
end module block_1lev_iterator

