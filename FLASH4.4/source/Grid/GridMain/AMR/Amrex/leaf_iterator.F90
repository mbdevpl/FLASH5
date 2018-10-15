!!****ih* source/Grid/GridMain/AMR/Amrex/leaf_iterator
!!
!! NAME
!!  leaf_iterator
!!
!! SYNOPSIS
!!
!!  use leaf_iterator, ONLY : leaf_iterator_t, build_iterator
!!  type(leaf_iterator_t) :: itor
!!  call build_iterator(leaf_iterator_t(OUT) :: itor,
!!                      integer(IN), optional :: level,
!!                      logical(IN), optional :: tiling )
!!
!! DESCRIPTION
!!
!! This Fortran module provides the leaf_iterator_t type, which can be
!! instantiated (or rather 'built') by iterator-using code outside of the
!! Grid unit to yield an iterator that iterates over leaf blocks only.
!!
!! Ideally, we will be able to use the AMReX iterator directly in the code
!! and client code will gain access to it through implementation-specific
!! code like Grid_getBlkIterator.
!!
!! This is a variant that uses type(block_1lev_iterator_t) as underlying iterator.
!!
!! ARGUMENTS
!!  itor     - the requested block iterator
!!  level    - iterate only over leaf blocks/tiles located at this level of
!!             refinement.
!!             A level value of UNSPEC_LEVEL is equivalent to omitting
!!             this optional argument.
!!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
!!             walk across all associated blocks on a tile-by-tile basis *if*
!!             the implementation supports this feature.  If a value is not
!!             given, is FALSE, or the implementation does not support tiling,
!!             the iterator will iterate on a block-by-block basis.
!!
!! NOTES
!!
!!  The notion of 'iteration over leaf blocks only' assumes that the
!!  AMReX Grid is organized in a valid octree structure. This does not
!!  mean that FLASH has to use the 'amrex_octree' module. I means,
!!  however, that the boxes that make up the multifabs that represent the
!!  FLASH solution variables like Unk, etc., satisfy certain rules; in
!!  particular, in a grid hierarchy the valid box of each block must be
!!  either fully covered or not covered at all by valid boxes of the next
!!  finer refinement level.
!!
!! SEE ALSO
!!  Grid_getLeafIterator
!!  Grid_releaseLeafIterator
!!  block_metadata_t
!!
!!****

#include "FortranLangFeatures.fh"
#include "constants.h"

module leaf_iterator

    use block_1lev_iterator, ONLY : block_1lev_iterator_t

    implicit none

    private

    public :: build_iterator, destroy_iterator

    !!****ic* leaf_iterator/leaf_iterator_t
    !!
    !! NAME
    !!  leaf_iterator_t
    !!
    !! DESCRIPTION
    !!  This class maintains a set of single-level iterators, which are used
    !!  internally to walk blocks/tiles.
    !!
    !!  NOTE: The three level integers as well as the index of li use FLASH's
    !!        1-based level indexing.
    !!****
    type, public :: leaf_iterator_t
        type(block_1lev_iterator_t), private, pointer :: li(:)
        integer                 :: first_level   = INVALID_LEVEL
        integer                 :: last_level    = INVALID_LEVEL
        integer                 :: level    = INVALID_LEVEL
        logical                 :: isValid = .FALSE.
        real,POINTER            :: fp(:,:,:,:)
    contains
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: blkMetaData
    end type leaf_iterator_t

    interface build_iterator
        procedure :: init_iterator
        procedure :: init_iterator_mfa
    end interface build_iterator

contains

    !!****im* leaf_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(leaf_iterator_t(OUT) :: itor,
    !!                 amrex_multifab(IN)    :: mfArray(:),
    !!                 integer(IN), optional :: level,
    !!                 logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the current AMReX octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!
    !! ARGUMENTS
    !!  itor     - the constructed iterator
    !!  mfArray  - an array of multfabs on which the iterator shall walk.  The
    !!             index is a 1-based index of the refinement levels.
    !!  level    - iterate only over LEAF blocks/tiles located at this level of
    !!             refinement.
    !!             A level value of UNSPEC_LEVEL is equivalent to omitting
    !!             this optional argument.
    !!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
    !!             walk across all associated blocks on a tile-by-tile basis *if*
    !!             the implementation supports this feature.  If a value is not
    !!             given, is FALSE, or the implementation does not support tiling,
    !!             the iterator will iterate on a block-by-block basis.
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine init_iterator_mfa(itor, mfArray, level, tiling)
      use amrex_multifab_module, ONLY : amrex_multifab

        type(leaf_iterator_t), intent(OUT) :: itor
        type(amrex_multifab),intent(IN),CONTIGUOUS :: mfArray(:)
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling

        integer :: lev, first, last
        logical :: v

        itor%level = UNSPEC_LEVEL
        if (present(level)) then
            first = level
            last = level
         else
            first = 1
            last = size(mfArray)
        end if

        allocate( itor%li (first : last) )

        itor%first_level = first
        itor%last_level = last
        itor%level = first

!!$        print*,'leaf_iterator_build: about to build 1lev iterators for this=',this%isValid,this%level,allocated(this%li)

        do lev=first,last
!!$           call amrex_mfiter_build(this%li(lev),mfArray(lev),tiling=tiling)
            itor%li(lev) = block_1lev_iterator_t(LEAF, mfArray(lev),lev,tiling=tiling)
!!$            call this%li( lev )%first()
            v = itor%li( lev )%is_valid()
            if (v .AND. .NOT. itor%isValid) then
               itor%isValid = .TRUE.
               itor%level   = lev
            end if
        end do

        if (.NOT. itor%isValid) then
           call destroy_iterator(itor)
        end if

!!$        print*,'leaf_iterator_build: done building 1lev iterators for this=',this%isValid,this%level,allocated(this%li)
!!$        call this%first()
      end subroutine init_iterator_mfa

    !!****im* leaf_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(leaf_iterator_t(OUT) :: itor,
    !!                 integer(IN), optional :: level,
    !!                 logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the current AMReX octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!
    !! ARGUMENTS
    !!  itor     - the constructed iterator
    !!  level    - iterate only over blocks/tiles located at this level of
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
    subroutine init_iterator(itor, level, tiling)
      use gr_physicalMultifabs,  ONLY : Unk

        type(leaf_iterator_t), intent(OUT)          :: itor
        integer,                intent(IN), optional :: level
        logical,                intent(IN), optional :: tiling

        integer :: lev, first, last
        integer :: finest_level
        logical :: v

        finest_level = size(Unk)

        itor%level = UNSPEC_LEVEL
        if (present(level)) then
            first = level
            last = level
         else
            first = 1
            last = finest_level
        end if

        allocate( itor%li (first : last) )

        itor%first_level = first
        itor%last_level = last
        itor%level = first


        do lev=first,last
            itor%li(lev) = block_1lev_iterator_t(LEAF,lev,tiling=tiling)
            v = itor%li( lev )%is_valid()
            if (v .AND. .NOT. itor%isValid) then
               itor%isValid = .TRUE.
               itor%level   = lev
            end if
        end do

        if (.NOT. itor%isValid) then
           call destroy_iterator(itor)
        end if
    end subroutine init_iterator

    !!****im* leaf_iterator_t/destroy_iterator
    !!
    !! NAME
    !!  destroy_iterator
    !!
    !! SYNPOSIS
    !!  Destroy given iterator
    !!
    !! DESCRIPTION
    !!  Clean-up block interator object at destruction
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(itor)
      type (leaf_iterator_t), intent(INOUT) :: itor

      integer :: lev

      if (associated(itor%li)) then
         do lev = itor%first_level, itor%last_level

            call itor%li(lev)%destroy_iterator()

         end do
         deallocate(itor%li)
      end if
      itor%isValid = .FALSE.

    end subroutine destroy_iterator

    !!****m* leaf_iterator_t/is_valid
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
        class(leaf_iterator_t), intent(IN) :: this
        logical :: ans

        ans = this%isValid
    end function is_valid

    !!****m* leaf_iterator_t/next
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
        class(leaf_iterator_t), intent(INOUT) :: this

        integer :: l
        logical :: v

!!$        print*,'leaf_iterator%next: about to do 1lev%next on this=',this%isValid,this%level,allocated(this%li)

        l = this%level

        call this%li( l )%next()
        v = this%li( l )%is_valid()

        do while (l .LT. this%last_level .AND. .NOT. v)
           l = l+1
!!$           call this%li( l )%first()   ! not necessary now!
           v = this%li( l )%is_valid()
        end do

        this%level = l
        this%isValid = v

    end subroutine next

    !!****m* leaf_iterator_t/blkMetaData
    !!
    !! NAME
    !!  blkMetaData
    !!
    !! SYNPOSIS
    !!  call itor%blkMetaData(block_metadata_t(OUT) : block)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block/tile currently set in the
    !!  iterator.
    !!
    !!****
    subroutine blkMetaData(this, blockDesc)
        use amrex_box_module,     ONLY : amrex_box

        use block_metadata,       ONLY : block_metadata_t

        class(leaf_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: blockDesc

        type(amrex_box) :: box, fabbox

        box = this%li( this%level )%tilebox()
        fabbox=this%li(this%level )%fabbox()

        blockDesc%grid_index        = this%li( this%level )%grid_index()
        blockDesc%level             = this%level

        ! FLASH uses 1-based spatial indices / AMReX uses 0-based
        blockDesc%limits(:, :) = 1
        blockDesc%limits(LOW,  1:N_DIM) = box%lo(1:N_DIM) + 1
        blockDesc%limits(HIGH, 1:N_DIM) = box%hi(1:N_DIM) + 1
        blockDesc%limitsGC(:, :) = 1
        blockDesc%limitsGC(LOW,  1:N_DIM) = fabbox%lo(1:N_DIM) + 1
        blockDesc%limitsGC(HIGH, 1:N_DIM) = fabbox%hi(1:N_DIM) + 1

        blockDesc%localLimits(LOW, :)   = blockDesc%limits(LOW, :)   - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimits(HIGH, :)  = blockDesc%limits(HIGH, :)  - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimitsGC(LOW, :) = 1
        blockDesc%localLimitsGC(HIGH, :)= blockDesc%limitsGC(HIGH, :)- blockDesc%limitsGC(LOW, :) + 1


    end subroutine blkMetaData

end module leaf_iterator

