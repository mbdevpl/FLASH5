!!****ih* source/Grid/GridMain/AMR/Amrex/block_1lev_iterator
!!
!! NAME
!!  block_1lev_iterator
!!
!! DESCRIPTION
!!  A class that defines an iterator facade around the AMReX MFIter
!!  (amrex_mfiter) such that client code may use the iterator for sequentially
!!  accessing specific blocks or tiles in the domain that exist only at a 
!!  single, given refinement level. 
!!
!!  Note that this iterator is meant only for internal FLASH use with 
!!  block_iterator_t.  No other code should need to use this code directly.
!!
!! SEE ALSO
!!  block_iterator_t
!!
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
    !! DESCRIPTION
    !!  For AMReX, this class is being used under the hood of the
    !!  block_iterator_t class and should not be used elsewhere.  As a result of
    !!  this restriction, for instance, there is no need to give access to
    !!  block_metadata_t structures.
    !!
    !!  Note that the value of level is specified using FLASH's 1-based level
    !!  indexing scheme.
    !!
    !!****
    type, public :: block_1lev_iterator_t
        type(amrex_mfiter),   private, pointer :: mfi      => NULL()
        type(amrex_multifab), private, pointer :: mf       => NULL()
        integer,              private          :: nodetype = LEAF
        integer,              private          :: level    = INVALID_LEVEL
        logical,              private          :: isValid  = .FALSE.
        integer,              private          :: finest_grid_level
        integer,              allocatable      :: dummy
    contains
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: grid_index
        procedure, public :: tilebox
        procedure, public :: fabbox
        procedure, public :: destroy_iterator
    end type block_1lev_iterator_t

    interface block_1lev_iterator_t
        procedure :: init_iterator
    end interface block_1lev_iterator_t

contains

    !!****im* block_1lev_iterator_t/block_1lev_iterator_t
    !!
    !! NAME
    !!  block_1lev_iterator_t
    !!
    !! SYNOPOSIS
    !!  itor = block_1lev_iterator_t(integer(IN)           :: nodetype,
    !!                               integer(IN)           :: level, 
    !!                               logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the the given refinement level.  The iterator is already
    !!  set to the first matching block/tile.
    !!
    !! ARGUMENTS
    !!  nodetype - the class of blocks to iterate over.  Acceptable values are
    !!             LEAF and ALL_BLKS.
    !!  level    - iterate only over blocks/tiles located at this level of
    !!             refinement.  Note that the level value must be given with
    !!             respect to FLASH's 1-based level index scheme.
    !!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
    !!             walk across all associated blocks on a tile-by-tile basis *if*
    !!             the implementation supports this feature.  If a value is not
    !!             given, is FALSE, or the implementation does not support tiling,
    !!             the iterator will iterate on a block-by-block basis.
    !!
    !! RETURN VALUE
    !!  The initialized iterator
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    function init_iterator(nodetype, level, tiling) result(this)
      use Driver_interface,      ONLY : Driver_abortFlash
      use gr_physicalMultifabs,  ONLY : unk
      use amrex_amrcore_module,  ONLY : amrex_get_finest_level

      integer, intent(IN)           :: nodetype
      integer, intent(IN)           :: level
      logical, intent(IN), optional :: tiling
      type(block_1lev_iterator_t)   :: this

      integer :: finest_level

      this%finest_grid_level = amrex_get_finest_level() ! 0-based finest existing level
      finest_level = this%finest_grid_level + 1 ! level and finest_level are 1-based
      if (level > finest_level) then
        call Driver_abortFlash("[init_iterator] No unk multifab for level")
      end if

      this%nodetype = nodetype
      this%level = level
      this%mf => unk(level-1)

      ! Don't permit repeated calls to init without intermediate destroy call
      if (associated(this%mfi)) then
        call Driver_abortFlash("[init_iterator] Destroy iterator before initializing again")
      end if

      allocate(this%mfi)
      call amrex_mfiter_build(this%mfi, this%mf, tiling=tiling)

      ! Set to True so that next() works
      this%isValid = .TRUE.

      ! Initial MFIter is not primed.  Advance to first compatible block.
      call this%next()
    end function init_iterator

    !!****im* block_1lev_iterator_t/destroy_iterator
    !!
    !! destroy_iterator
    !!
    !! SYNPOSIS
    !!  itor%destroy_iterator
    !!
    !! DESCRIPTION
    !!  Clean-up block iterator and internal resources.  Note that this is not a
    !!  destructor and therefore must be called manually.
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(this)
      class(block_1lev_iterator_t), intent(INOUT) :: this

      if (associated(this%mfi)) then
        call amrex_mfiter_destroy(this%mfi)
        deallocate(this%mfi)
        nullify(this%mfi)
      end if

      nullify(this%mf)
      this%isValid = .FALSE.
    end subroutine destroy_iterator

    !!****m* block_1lev_iterator_t/is_valid
    !!
    !! NAME
    !!  is_valid
    !!
    !! SYNPOSIS
    !!  logical valid = itor%is_valid()
    !!
    !! DESCRIPTION
    !!  Determine if the iterator is currently set to a valid block/tile.
    !!
    !! RETURN VALUE 
    !!  True if iterator is currently set to a valid block/tile
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
    !!  Advance the iterator to the next block/tile managed by process and
    !!  that meets the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next(this)
        use amrex_box_module, ONLY : amrex_box
        use Driver_interface, ONLY : Driver_abortFlash

        class(block_1lev_iterator_t), intent(INOUT) :: this

        type(amrex_box) :: bx
        logical         :: hasChildren

        if (this%isValid) then
           do
              this%isValid = this%mfi%next()
              if (.NOT. this%isValid) then
                 exit
              else
                 select case (this%nodetype)
                 case(ALL_BLKS)
                    exit
                 case(LEAF)
                    bx = this%mfi%tilebox()
                    hasChildren = boxIsCovered(bx, this%level-1, this%finest_grid_level)
                    if (.NOT.hasChildren) exit
                 case default
                    call Driver_abortFlash("[block_1lev_iterator]: Unsupported nodetype")
                 end select
              end if
           end do
        else
           call Driver_abortFlash("[block_1lev_iterator]: attempting next() on invalid!")
        end if

    contains

        logical function boxIsCovered(bx,lev,finest_level) result(covered)
          use amrex_boxarray_module, ONLY : amrex_boxarray
          use amrex_amrcore_module,  ONLY : amrex_ref_ratio, &
                                            amrex_get_boxarray

          !IMPORTANT: data in bx is changed on return!
          type(amrex_box), intent(INOUT) :: bx
          integer,         intent(IN)    :: lev
          integer,         intent(IN)    :: finest_level ! Passing this saves a function call.

          type(amrex_boxarray) :: fba
          integer :: rr

          ! Assume lev is 0-based
          if (lev .GE. finest_level) then
             covered = .FALSE.
          else
             fba = amrex_get_boxarray(lev+1)
             rr = amrex_ref_ratio(lev)

             ! Note: this modifies bx, do not use naively after this!
             call bx%refine(rr)
             covered = fba%intersects(bx)
          end if
        end function boxIsCovered

    end subroutine next

    !!****m* block_1lev_iterator_t/grid_index
    !!
    !! NAME
    !!  grid_index
    !!
    !! SYNPOSIS
    !!  idx = itor%grid_index()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next block managed by process and that meets
    !!  the iterator constraints given at instantiation.
    !!
    !!****
    function grid_index(this) result(idx)
      class(block_1lev_iterator_t), intent(IN) :: this
      integer                                  :: idx

      idx = this%mfi%grid_index()
    end function grid_index

    !!****m* block_1lev_iterator_t/tilebox
    !!
    !! NAME
    !!  tilebox
    !!
    !! SYNPOSIS
    !!  box = itor%tilebox()
    !!
    !! DESCRIPTION
    !!  Obtain the box without guardcells of the block/tile currently
    !!  "loaded" into the iterator.
    !!
    !! RETURN VALUE
    !!  An AMReX box object.  The index space of the box is the index
    !!  space of the multifab used to construct the underlying MFIter.
    !!  The spatial indices of the box use AMReX's 0-based scheme.
    !!
    !!****
    function tilebox(this) result(bx)
      use amrex_box_module, ONLY : amrex_box

      class(block_1lev_iterator_t), intent(in) :: this
      type(amrex_box)                          :: bx

      bx = this%mfi%tilebox()
    end function tilebox

    !!****m* block_1lev_iterator_t/fabbox
    !!
    !! NAME
    !!  fabbox
    !!
    !! SYNPOSIS
    !!  box = itor%fabbox()
    !!
    !! DESCRIPTION
    !!  Obtain the box wth guardcells of the block/tile currently
    !!  "loaded" into the iterator.
    !!
    !! RETURN VALUE
    !!  An AMReX box object.  The index space of the box is the index
    !!  space of the multifab used to construct the underlying MFIter.
    !!  The spatial indices of the box use AMReX's 0-based scheme.
    !!
    !!****
    function fabbox(this) result(bx)
      use amrex_box_module, ONLY : amrex_box

      class(block_1lev_iterator_t), intent(in) :: this
      type(amrex_box)                          :: bx
      
      bx = this%mfi%fabbox()
    end function fabbox

end module block_1lev_iterator

