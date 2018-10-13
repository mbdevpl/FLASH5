!!****ih* source/Grid/GridMain/AMR/Amrex/gr_iterator
!!
!! NAME
!!  gr_iterator
!!
!! DESCRIPTION
!!  A class that defines a full-featured iterator for sequentially accessing
!!  specific blocks or tiles in the physical domain.  At initialization, the
!!  client code informs the initialization routine what blocks/tiles need to
!!  be accessed by the client code via the iterator.
!!
!!  Please refer to the documentation of gr_getBlkIterator for more
!!  information regarding iterator initialization.
!!
!!  NOTE: This iterator is meant for use only within the Grid unit.
!!
!! EXAMPLE
!!  The following example demonstrates looping over all blocks defined on
!!  the coarsest level.
!!
!!  type(block_metadata_t) :: blockDesc 
!!
!!  call gr_getBlkIterator(itor, level=1)
!!  do while (itor%is_valid())
!!    call itor%blkMetaData(blockDesc)
!!
!!    call Grid_getBlkPtr(blockDesc, solnData, CENTER)
!!                          ...
!!      work with cell-centered data of current block
!!                          ... 
!!    call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)
!!
!!    call itor%next()
!!  end do
!!  call gr_releaseBlkIterator(itor)
!!
!! SEE ALSO
!!  gr_getBlkIterator
!!  gr_releaseBlkIterator
!!  block_descriptor_t
!!
!!****

#include "FortranLangFeatures.fh"
#include "constants.h"
#include "Flash.h"

module gr_iterator

    use block_1lev_iterator, ONLY : block_1lev_iterator_t

    implicit none

    private

    public :: build_iterator, destroy_iterator

    !!****ic* gr_iterator/gr_iterator_t
    !!
    !! NAME
    !!  gr_iterator_t
    !!
    !! DESCRIPTION
    !!  This class maintains a set of single-level iterators, which are used
    !!  internally to walk blocks/tiles.
    !!
    !!  NOTE: The three level integers as well as the index of li use FLASH's
    !!        1-based level indexing.
    !!****
    type, public :: gr_iterator_t
        type(block_1lev_iterator_t), private, pointer :: li(:)
        integer,                     private              :: first_level = INVALID_LEVEL
        integer,                     private              :: last_level  = INVALID_LEVEL
        integer,                     private              :: level       = INVALID_LEVEL
        logical,                     private              :: isValid     = .FALSE.
    contains
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: blkMetaData
    end type gr_iterator_t

    interface build_iterator
        procedure :: init_iterator
    end interface build_iterator

contains

    !!****im* gr_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(gr_iterator_t(OUT)    :: itor,
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
    !!        interface --- gr_getBlkIterator/gr_releaseBlkIterator.
    !!
    !! ARGUMENTS
    !!  itor     - the constructed iterator
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS).
    !!             Refer to the documentation for the AMReX version of
    !!             gr_getBlkIterator for more information.
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
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine init_iterator(itor, nodetype, level, tiling)
      use amrex_amrcore_module,  ONLY : amrex_get_finest_level

      type(gr_iterator_t), intent(OUT)          :: itor
      integer,             intent(IN)           :: nodetype
      integer,             intent(IN), optional :: level
      logical,             intent(IN), optional :: tiling

      integer :: lev
      integer :: finest_level
      logical :: is_lev_valid
            
      finest_level = amrex_get_finest_level() + 1

      associate(first => itor%first_level, &
                last  => itor%last_level)
        if (present(level)) then
           if (level .NE. UNSPEC_LEVEL) then
            ! Construct do nothing iterator if no blocks on level
              if (level > finest_level) then
                 itor%isValid = .FALSE.
                 RETURN
              end if

              first = level
              last = level
           else
              first = 1
              last = finest_level
           end if
        else
            first = 1
            last = finest_level
        end if
        itor%level = first
 
        allocate( itor%li(first : last) )

        do lev=first, last
            itor%li(lev) = block_1lev_iterator_t(nodetype, lev, tiling=tiling)
            is_lev_valid = itor%li(lev)%is_valid()
            if (is_lev_valid .AND. .NOT. itor%isValid) then
               itor%isValid = .TRUE.
               itor%level   = lev
            end if
        end do
      end associate
    end subroutine init_iterator

    !!****im* gr_iterator_t/destroy_iterator
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
      type (gr_iterator_t), intent(INOUT) :: itor

      integer :: lev

      if (associated(itor%li)) then
         do lev = itor%first_level, itor%last_level

            call itor%li(lev)%destroy_iterator()

         end do
         deallocate(itor%li)
      end if
      itor%isValid = .FALSE.
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
    function is_valid(this) result(ans)
        class(gr_iterator_t), intent(IN) :: this
        logical :: ans

        ans = this%isValid
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

            this%isValid = is_li_valid
        end associate
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
    !!  Obtain meta data that characterizes the block/tile currently set in the
    !!  iterator.
    !!
    !!****
    subroutine blkMetaData(this, blockDesc)
        use amrex_box_module, ONLY : amrex_box
        use block_metadata,   ONLY : block_metadata_t

        class(gr_iterator_t),   intent(IN)  :: this
        type(block_metadata_t), intent(OUT) :: blockDesc

        type(amrex_box) :: box, fabbox
       
        box    = this%li( this%level )%tilebox()
        fabbox = this%li( this%level )%fabbox()

        blockDesc%grid_index = this%li( this%level )%grid_index()
        blockDesc%level      = this%level
        
        ! FLASH uses 1-based spatial indices / AMReX uses 0-based
        blockDesc%limits(:, :) = 1
        blockDesc%limits(LOW,  1:NDIM) = box%lo(1:NDIM) + 1
        blockDesc%limits(HIGH, 1:NDIM) = box%hi(1:NDIM) + 1
        blockDesc%limitsGC(:, :) = 1
        blockDesc%limitsGC(LOW,  1:NDIM) = fabbox%lo(1:NDIM) + 1
        blockDesc%limitsGC(HIGH, 1:NDIM) = fabbox%hi(1:NDIM) + 1

        blockDesc%localLimits(LOW, :)    =   blockDesc%limits(LOW, :) &
                                           - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimits(HIGH, :)   =   blockDesc%limits(HIGH, :) &
                                           - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimitsGC(LOW, :)  = 1
        blockDesc%localLimitsGC(HIGH, :) =   blockDesc%limitsGC(HIGH, :) &
                                           - blockDesc%limitsGC(LOW, :) + 1
    end subroutine blkMetaData
 
end module gr_iterator

