!!****ih* source/Grid/GridMain/AMR/paramesh/block_iterator
!!
!!
!!
!!****

module block_iterator

#include "Flash.h"
    use tree, ONLY : lnblocks, lrefine, lrefine_max
    use Grid_interface, ONLY : Grid_getBlkPtr
    use Grid_interface, ONLY : Grid_getBlkIndexLimits
    use Grid_interface, ONLY : Grid_getBlkCornerID
    use Grid_interface, ONLY : Grid_blockMatch

    implicit none

#define IMPURE_ELEMENTAL
#define CONTIGUOUS_POINTER pointer
#include "constants.h"
    private
  
    integer,parameter :: ndims=N_DIM
    integer,parameter :: amrex_real=kind(1.0)

    !!****ic* block_iterator/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !!****
    type, public :: block_iterator_t
        integer :: cur         = 1
        integer :: nodetype    = LEAF 
        integer :: gds         = CENTER
        integer :: lev         = INVALID_LEVEL
        integer :: cellIdxBase = -2
        ! different variants for cell index numbering:
        !  1 = normal FLASH convention:     per block, leftmost guard cell = 1
        !  0 = zero-based FLASH convention: per block, leftmost guard cell = 0
        ! -1 = global convention for a refinement level: leftmost guard cell = 1
        ! -2 = global convention for a refinement level: leftmost inner cell = 1
    contains
        procedure, public :: clear
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: blkLimits
        procedure, public :: blkID
        procedure, public :: blkCornerID
        procedure, public :: blkStride
        procedure, public :: blkLevel
        procedure, public :: blkDataPtr
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
    !!                                           integer(IN)         :: gds,
    !!                                           level(IN), optional :: level)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks
    !!  within the current paramesh octree structure.  The iterator is already
    !!  set to the first matching block.
    !!
    !! ARGUMENTS
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS)
    !!  gds      - integer value specifying data structure.  The options are 
    !!             defined in constants.h (e.g. CENTER, FACEX, SCRATCH_CTR)
    !!  level    - if nodetype is LEAF, PARENT, ANCESTOR, or REFINEMENT, then 
    !!             iterate only over blocks located at this level of 
    !!             octree structure
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    function init_iterator(nodetype, gds, level) result(this)
        integer, intent(IN)           :: nodetype
        integer, intent(IN)           :: gds 
        integer, intent(IN), optional :: level
        type(block_iterator_t)        :: this
 
        this%nodetype = nodetype
        this%gds = gds
        this%lev = INVALID_LEVEL 
        if (present(level)) then
            this%lev = level
        end if

        ! Search for the first valid block
        this%cur = 0
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

        call this%clear()
    end subroutine destroy_iterator
#endif

    !!****m* block_iterator_t/clear
    !!
    !! NAME
    !!  clear 
    !!
    !! SYNPOSIS
    !!  call itor%clear() 
    !!
    !! DESCRIPTION
    !!  Reset iterator to the initial block managed by process
    !!
    !!****
    subroutine clear(this)
        class(block_iterator_t), intent(INOUT) :: this

        this%cur = 0
        call this%next()
    end subroutine clear
 
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
        use tree, ONLY : lnblocks

        class(block_iterator_t), intent(IN) :: this

        is_valid = (this%cur <= lnblocks)
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

        integer :: j = 0
  
        ! DEVNOTE: Move this check inside of Grid_blockMatch
        if (this%lev == INVALID_LEVEL) then
            do j = this%cur + 1, lnblocks
                if (Grid_blockMatch(j, this%nodetype)) EXIT
            end do
        else
            do j = this%cur + 1, lnblocks
                if (Grid_blockMatch(j, this%nodetype, this%lev)) EXIT
            end do
        end if

        this%cur = j
    end subroutine next

    !!****m* block_iterator_t/blkLimits
    !!
    !! NAME
    !!  blkLimits
    !!
    !! SYNPOSIS
    !!  call itor%blkLimits(integer(INOUT) :: limits(LOW:HIGH, MDIM))
    !!
    !! DESCRIPTION
    !!  DEVNOTE: WRITE THIS AND EXPLAIN DIFFERENCES BY cellIdxBase
    !!
    !! ARGUMENTS
    !!  limits - the limits
    !!
    !! SEE ALSO
    !!  Grid_getBlkIndexLimits
    !!  Grid_getBlkCornerID
    !!
    !!****
    subroutine blkLimits(this, limits)
        use tree, ONLY : lrefine, lrefine_max

        class(block_iterator_t), intent(IN)    :: this
        integer,                 intent(INOUT) :: limits(LOW:HIGH, MDIM)
        
        integer, dimension(MDIM)           :: cornerID, strideUnused
        integer, dimension(LOW:HIGH, MDIM) :: blkLim, blkLimGC

        associate(lo => limits(LOW, :), hi => limits(HIGH, :), blkId => this%cur)
            call Grid_getBlkIndexLimits(blkID, blkLim, blkLimGC)
            lo = blkLim(LOW ,:)                   !DEV: or blkLimGC(LOW ,:) ?
            hi = blkLim(HIGH,:)                   !DEV: or blkLimGC(HIGH,:) ?
            if (this%cellIdxBase==-1) then
               call Grid_getBlkCornerID(blkID, cornerID, strideUnused)
               cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(blkID)) + 1
               lo(:) = lo(:) - 1 + cornerID(:)
               hi(:) = hi(:) - 1 + cornerID(:)
            else if (this%cellIdxBase==-2) then
               call Grid_getBlkCornerID(blkID, cornerID, strideUnused)
               cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(blkID)) + 1
               lo(:) = lo(:) - 1 + cornerID(:)
               hi(:) = hi(:) - 1 + cornerID(:)
               lo(1:ndims) = lo(1:ndims) - NGUARD
               hi(1:ndims) = hi(1:ndims) - NGUARD
            else if (this%cellIdxBase==0) then
               lo(:) = lo(:) - 1
               hi(:) = hi(:) - 1
            end if
        end associate
    end subroutine blkLimits

    !!****m* block_iterator_t/blkID
    !!
    !! NAME
    !!  blkID
    !!
    !! SYNPOSIS
    !!  integer id = itor%blkID()
    !!
    !! DESCRIPTION
    !!  Obtain the paramesh-based integer that uniquely indexes the block
    !!  currently set in the iterator.
    !!
    !! RETURN VALUES
    !!  The integer index for the block
    !!
    !!****
    integer function blkID(this)
        class(block_iterator_t), intent(IN) :: this

        blkID = this%cur
    end function blkID

    function blkCornerID(this) result(cid)
        class(block_iterator_t), intent(IN) :: this
        ! TODO: What is dimension of this
        integer :: cid(MDIM)

        cid = 0
    end function blkCornerID

    function blkStride(this) result(cid)
        class(block_iterator_t), intent(IN) :: this
        ! TODO: What is dimension of this
        integer :: cid(MDIM)
        
        cid = 0
    end function blkStride

    !!****m* block_iterator_t/blkLevel
    !!
    !! NAME
    !!  blkLevel
    !!
    !! SYNPOSIS
    !!  integer level = itor%blkLevel()
    !!
    !! DESCRIPTION
    !!  Obtain the octree level of the block currently set in the iterator.
    !!
    !! RETURN VALUE
    !!  The octree level 
    !!
    !!****
    integer function blkLevel(this)
        use tree, ONLY : lrefine

        class(block_iterator_t), intent(IN) :: this

        blkLevel = lrefine(this%cur)
    end function blkLevel 

    !!****m* block_iterator_t/blkDataPtr 
    !!
    !! NAME
    !!  blkDataPtr
    !!
    !! SYNPOSIS
    !!  real ptr(:,:,:,:) = itor%blkDataPtr()
    !!
    !! DESCRIPTION
    !!  Obtain a pointer to the full set of data defined on the block currently 
    !!  set in the iterator.  The extent of the data is dependent on the gds
    !!  value given to the iterator at instantiation.
    !!
    !! RETURN VALUE
    !!  The pointer
    !!
    !! NOTES
    !!  The array pointed to is subject to having its index order reordered.
    !!
    !!****
    function blkDataPtr(this) result(ptr)
        use tree, ONLY : lrefine, lrefine_max

        class(block_iterator_t), intent(IN)                      :: this
        real(amrex_real), CONTIGUOUS_POINTER, dimension(:,:,:,:) :: ptr
        
        real(amrex_real), CONTIGUOUS_POINTER, dimension(:,:,:,:) :: fp
        integer, dimension(LOW:HIGH,MDIM)    :: blkLim ,blkLimGC
        integer, dimension(MDIM)             :: cornerID, strideUnused
        integer, dimension(MDIM)             :: lo, hi
        
        associate (blkID => this%cur)
            call Grid_getBlkIndexLimits(blkID, blkLim, blkLimGC)
            call Grid_getBlkPtr(blkID, fp, this%gds)

            lo = blkLimGC(LOW ,:)                 !DEV: or blkLimGC(LOW ,:) ?
            if (this%cellIdxBase==-1) then
               call Grid_getBlkCornerID(blkID, cornerID, strideUnused)
               cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(blkID)) + 1
               lo(:) = lo(:) - 1 + cornerID(:)
            else if (this%cellIdxBase==-2) then
               call Grid_getBlkCornerID(blkID,cornerID,strideUnused)
               cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(blkID)) + 1
               lo(:) = lo(:) - 1 + cornerID(:)
               lo(1:ndims) = lo(1:ndims) - NGUARD
            else if (this%cellIdxBase==0) then
               lo(:) = lo(:) - 1
            end if
#ifdef INDEXREORDER
            ptr(lo(1):, lo(2):, lo(3):, 1:) => fp
#else
            ptr(1:, lo(1):, lo(2):, lo(3):) => fp
#endif
        end associate
    end function blkDataPtr

end module block_iterator

