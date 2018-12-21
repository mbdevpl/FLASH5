!!****ih* source/Grid/GridMain/AMR/Amrex/flash_iterator
!!
!! NAME
!!  flash_iterator
!!
!! DESCRIPTION
!!  A class that defines a full-featured iterator for sequentially accessing
!!  specific blocks or tiles in the physical domain.  At initialization, the
!!  client code informs the initialization routine what blocks/tiles need to
!!  be accessed by the client code via the iterator.
!!
!!  Please refer to the documentation of Grid_getTileIterator for more
!!  information regarding iterator initialization.
!!
!! EXAMPLE
!!  The following example demonstrates looping over all blocks defined on
!!  the coarsest level.
!!
!!  type(flash_tile_t) :: tileDesc
!!
!!  call Grid_getTileIterator(itor, ALL_BLKS, level=1)
!!  do while (itor%isValid())
!!    call itor%currentTile(tileDesc)
!!
!!    call tileDesc%getDataPtr(solnData, CENTER)
!!                          ...
!!      work with cell-centered data of current block
!!                          ... 
!!    call tileDesc%releaseDataPtr(solnData, CENTER)
!!
!!    call itor%next()
!!  end do
!!  call Grid_releaseTileIterator(itor)
!!
!! SEE ALSO
!!  Grid_getTileIterator
!!  Grid_releaseTileIterator
!!  flash_tile_t
!!
!!****

#include "FortranLangFeatures.fh"
#include "constants.h"
#include "Flash.h"

module flash_iterator
    use block_1lev_iterator, ONLY : block_1lev_iterator_t

    implicit none

    private

    public :: build_iterator, destroy_iterator

    !!****ic* flash_iterator/flash_iterator_t
    !!
    !! NAME
    !!  flash_iterator_t
    !!
    !! DESCRIPTION
    !!  This class maintains a set of single-level iterators, which are used
    !!  internally to walk blocks/tiles.
    !!
    !!  NOTE: The three level integers as well as the index of li use FLASH's
    !!        1-based level indexing.
    !!****
    type, public :: flash_iterator_t
        type(block_1lev_iterator_t), private, pointer :: li(:)       => null()
        integer,                     private          :: first_level = INVALID_LEVEL
        integer,                     private          :: last_level  = INVALID_LEVEL
        integer,                     private          :: level       = INVALID_LEVEL
        logical,                     private          :: is_valid    = .FALSE.
    contains
        procedure, public :: isValid
        procedure, public :: next
        procedure, public :: currentTile
    end type flash_iterator_t

    interface build_iterator
        procedure :: init_iterator
    end interface build_iterator

contains

    !!****im* flash_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(flash_iterator_t(OUT) :: itor,
    !!                 integer(IN)           :: nodetype,
    !!                 integer(IN), optional :: level, 
    !!                 logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the current AMReX octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!
    !!  NOTE: Prefer iterator acquisition/destruction via Grid unit local 
    !!        interface --- Grid_getTileIterator/Grid_releaseTileIterator.
    !!
    !! ARGUMENTS
    !!  itor     - the constructed iterator
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS).
    !!             Refer to the documentation for the AMReX version of
    !!             Grid_getTileIterator for more information.
    !!  level    - iterate only over all blocks/tiles of the correct nodetype
    !!             that are located at this level of refinement.  Note that the
    !!             level value must be given with respect to FLASH's 1-based
    !!             level index scheme.  If no level value is given, then
    !!             iteration is not restricted to any level.
    !!             A level value of UNSPEC_LEVEL is equivalent to omitting
    !!             this optional argument.
    !!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
    !!             walk across all associated blocks on a tile-by-tile basis *if*
    !!             the implementation supports this feature.  If a value is not
    !!             given, is FALSE, or the implementation does not support tiling,
    !!             the iterator will iterate on a block-by-block basis.
    !!
    !!****
    subroutine init_iterator(itor, nodetype, level, tiling, tileSize)
      use amrex_amrcore_module, ONLY : amrex_get_finest_level

      type(flash_iterator_t), intent(OUT) :: itor
      integer,                intent(IN)  :: nodetype
      integer,                intent(IN)  :: level
      logical,                intent(IN)  :: tiling
      integer,                intent(IN)  :: tileSize(1:MDIM)

      integer :: lev
      integer :: finest_level
      logical :: is_lev_valid

      finest_level = amrex_get_finest_level() + 1

      associate(first => itor%first_level, &
                last  => itor%last_level)
        if (level .NE. UNSPEC_LEVEL) then
         ! Construct do nothing iterator if no blocks on level
           if (level > finest_level) then
              itor%is_valid = .FALSE.
              RETURN
           end if

           first = level
           last = level
        else
           first = 1
           last = finest_level
        end if
        itor%level = first
 
        allocate( itor%li(first : last) )

        do lev=first, last
            itor%li(lev) = block_1lev_iterator_t(nodetype, lev, &
                                                 tiling=tiling, &
                                                 tileSize=tileSize)
            is_lev_valid = itor%li(lev)%is_valid()
            if (is_lev_valid .AND. .NOT. itor%is_valid) then
               itor%is_valid = .TRUE.
               itor%level   = lev
            end if
        end do
      end associate
    end subroutine init_iterator

    !!****im* flash_iterator_t/destroy_iterator
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
      type (flash_iterator_t), intent(INOUT) :: itor

      integer :: lev

      if (associated(itor%li)) then
         do lev = itor%first_level, itor%last_level
            call itor%li(lev)%destroy_iterator()
         end do

         deallocate(itor%li)
         nullify(itor%li)
      end if

      itor%is_valid = .FALSE.
    end subroutine destroy_iterator

    !!****m* flash_iterator_t/isValid
    !!
    !! NAME
    !!  isValid
    !!
    !! SYNPOSIS
    !!  logical valid = itor%isValid()
    !!
    !! DESCRIPTION
    !!  Determine if the iterator is currently set to a valid block/tile.
    !!
    !! RETURN VALUE 
    !!  True if iterator is currently set to a valid block/tile.
    !!
    !!****
    function isValid(this)
        class(flash_iterator_t), intent(IN) :: this
        logical :: isValid

        isValid = this%is_valid
    end function isValid

    !!****m* flash_iterator_t/next
    !!
    !! NAME
    !!  next
    !!
    !! SYNPOSIS
    !!  call itor%next()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next block/tile managed by process and that meets
    !!  the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next(this)
        class(flash_iterator_t), intent(INOUT) :: this

        logical :: is_li_valid

        associate(lev => this%level)
            call this%li( lev )%next()
            is_li_valid = this%li( lev )%is_valid()

            ! Search for next allowable level that has blocks meeting our
            ! criteria
            do while ((lev .LT. this%last_level) .AND. (.NOT. is_li_valid))
               lev = lev + 1
               is_li_valid = this%li( lev )%is_valid()
            end do

            this%is_valid = is_li_valid
        end associate
    end subroutine next

    !!****m* flash_iterator_t/currentTile
    !!
    !! NAME
    !!  currentTile 
    !!
    !! SYNPOSIS
    !!  call itor%currentTile(flash_tile_t(OUT) : block)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block/tile currently set in the
    !!  iterator.
    !!
    !!****
    subroutine currentTile(this, tileDesc)
        use amrex_box_module, ONLY : amrex_box
        use flash_tile,       ONLY : flash_tile_t

        class(flash_iterator_t), intent(IN)  :: this
        type(flash_tile_t),      intent(OUT) :: tileDesc

        type(amrex_box) :: box

        tileDesc%level      = this%level
        tileDesc%grid_index = this%li( this%level )%grid_index()
        tileDesc%tile_index = this%li( this%level )%local_tile_index()

        ! FLASH uses 1-based spatial indices / AMReX uses 0-based
        box = this%li( this%level )%tilebox()
        tileDesc%limits(:, :) = 1
        tileDesc%limits(LOW,  1:NDIM) = box%lo(1:NDIM) + 1
        tileDesc%limits(HIGH, 1:NDIM) = box%hi(1:NDIM) + 1

        box = this%li( this%level )%growntilebox()
        tileDesc%limitsGC(:, :) = 1
        tileDesc%limitsGC(LOW,  1:NDIM) = box%lo(1:NDIM) + 1
        tileDesc%limitsGC(HIGH, 1:NDIM) = box%hi(1:NDIM) + 1

        box = this%li( this%level )%fabbox()
        tileDesc%blkLimitsGC(:, :) = 1
        tileDesc%blkLimitsGC(LOW,  1:NDIM) = box%lo(1:NDIM) + 1
        tileDesc%blkLimitsGC(HIGH, 1:NDIM) = box%hi(1:NDIM) + 1
    end subroutine currentTile 
 
end module flash_iterator

