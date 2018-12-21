!!****ih* source/Grid/GridMain/AMR/Amrex/flash_tile
!!
!!  NOTES
!!    It is generally expected that these objects will never be created directly
!!    in FLASH and therefore there are no "constructors" specified.  In the
!!    even that low-level code in the Grid unit should need to create an object
!!    it can be done by directly setting all public fields to appropriate
!!    values.
!!
!!    DEV: TODO Fix this documentation.
!!    The data contained in the fields of this derived type are
!!
!!    level - the refinement level at which the block resides.  Levels are
!!            referenced with the index set {1, 2, ..., gr_lRefineMax}, where
!!            1 refers to the coarsest level and mesh refinement increases with
!!            increasing level index.  The value of gr_lRefineMax is fixed to
!!            the value of the runtime parameter lrefine_max.
!!
!!    grid_index - used internally by AMReX.  Client code should never make
!!                 direct use of the value stored in this field.
!!
!!    tile_index - used internally by AMReX.  Client code should never make
!!                 direct use of the value stored in this field.
!!
!!    limits - a two-dimensional array where limits(LOW, :) is the array
!!             containing the (i, j, k) coordinate in index space of the 
!!             lower-leftmost corner cell of the block.  Similarly, 
!!             limits(HIGH, :) is the array containing the (i, j, k) 
!!             coordinate of the upper-rightmost corner cell of the block.
!!             The two coordinates are specified for the cell-center index space
!!             defined so that 
!!                - i is an integer in {1, 2, ..., n_x_cells},
!!                - j is an integer in {1, 2, ..., n_y_cells, and
!!                - k is an integer in {1, 2, ..., n_z_cells},
!!             where for FLASH, 
!!                  n_[xyz]_cells = 2^(level - 1) * N[XYZ]B * nblock[xyz].
!!             This index space is the global index space so that (1, 1, 1) is
!!             the lower-leftmost cell in the full domain of the problem.
!!
!!    blkLimitsGC - a two-dimensional array identical to limits except that the
!!                  two cell coordinates in index space refer to the lower-leftmost
!!                  and the upper-rightmost guardcells for the tile's subblock.
!!                  Note that if the tile is in fact a block, then blkLimitsGC 
!!                  has the same contents as limitsGC.  The coordinates are
!!                  specified w.r.t. an expanded version of the global
!!                  cell-center index space used for limits so that
!!                   - i is an integer in {1-NGUARD, ..., n_x_cells+NGUARD},
!!                   - j is an integer in {1-NGUARD, ..., n_y_cells+NGUARD},
!!                   - k is an integer in {1-NGUARD, ..., n_z_cells+NGUARD}.
!!                  Hence, (1, 1, 1) still refers to the lower-leftmost guardcell for
!!                  the full domain of the problem.  If gr_[ijk]guard=4, then 
!!                  (-3, -3, -3) indexes the lower-leftmost guardcell in the
!!                  full domain.
!!
!!    limitsGC - a two-dimensional array identical to limits except that the
!!               two cell coordinates in index space refer to the lower-leftmost
!!               and the upper-rightmost guardcells for the block.  The
!!               coordinates are specified w.r.t. an expanded version of the
!!               global cell-center index space used for limits so that
!!                - i is an integer in {1-NGUARD, ..., n_x_cells+NGUARD},
!!                - j is an integer in {1-NGUARD, ..., n_y_cells+NGUARD},
!!                - k is an integer in {1-NGUARD, ..., n_z_cells+NGUARD}.
!!               Hence, (1, 1, 1) still refers to the lower-leftmost guardcell for
!!               the full domain of the problem.  If gr_[ijk]guard=4, then 
!!               (-3, -3, -3) indexes the lower-leftmost guardcell in the
!!               full domain.
!!
!!****

#include "constants.h"
#include "Flash.h"

module flash_tile
    implicit none

    private

    type, public :: flash_tile_t
        integer, public :: level
        integer, public :: grid_index
        integer, public :: tile_index
        integer, public :: limits(LOW:HIGH, MDIM)
        integer, public :: limitsGC(LOW:HIGH, MDIM)
        integer, public :: blkLimitsGC(LOW:HIGH, MDIM)
    contains
        procedure, public :: deltas
        procedure, public :: boundBox
        procedure, public :: coordinates
!        procedure, public :: faceAreas
!        procedure, public :: cellVolumes
!        procedure, public :: physicalSize
!        procedure, public :: onDomainBoundary
!        procedure, public :: outsideDomain
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
        procedure, public :: enclosingBlock
    end type flash_tile_t

contains

    subroutine deltas(this, dx)
        use amrex_amrcore_module, ONLY : amrex_geom

        class(flash_tile_t), intent(IN)  :: this
        real,                intent(OUT) :: dx(1:MDIM)

        ! AMReX uses zero-based level indexing, but FLASH assumes one-based
        dx(:) = 0.0
        dx(1:MDIM) = amrex_geom(this%level - 1)%dx(1:MDIM)
    end subroutine deltas

    subroutine boundBox(this, box)
        use amrex_amrcore_module,  ONLY : amrex_geom
        use amrex_geometry_module, ONLY : amrex_problo

        class(flash_tile_t), intent(IN)  :: this
        real,                intent(OUT) :: box(LOW:HIGH, 1:MDIM)

        ! DEV: FIXME How to manage matching amrex_real to FLASH real
        box(:, :) = 1.0
        associate(x0 => amrex_problo, &
                  dx => amrex_geom(this%level - 1)%dx, &
                  lo => this%limits(LOW,  :), &
                  hi => this%limits(HIGH, :))
!            ! lo is 1-based cell-index of lower-left cell in block 
!            ! hi is 1-based cell-index of upper-right cell in block
            box(LOW,  1:NDIM) = x0(1:NDIM) + (lo(1:NDIM) - 1)*dx(1:NDIM)
            box(HIGH, 1:NDIM) = x0(1:NDIM) + (hi(1:NDIM)    )*dx(1:NDIM)
        end associate
    end subroutine boundBox 

    subroutine coordinates(this, axis, edge, guardcell, coords)
        use amrex_amrcore_module,  ONLY : amrex_geom
        use amrex_geometry_module, ONLY : amrex_problo

        use Driver_interface,      ONLY : Driver_abortFlash

        class(flash_tile_t), intent(IN)  :: this
        integer,             intent(IN)  :: axis
        integer,             intent(IN)  :: edge
        logical,             intent(IN)  :: guardcell
        real,                intent(OUT) :: coords(:)

        real    :: shift
        integer :: i

#ifdef DEBUG_GRID
        integer :: width

        if((axis/=IAXIS) .and. (axis/=JAXIS) .and. (axis/=KAXIS)) then
           call Driver_abortFlash("Get Coords : invalid axis, must be IAXIS, JAXIS or KAXIS ")
        end if

        width = this%limits(HIGH, axis) - this%limits(LOW, axis) + 1
        if (guardcell) then
            width = width + 2*NGUARD
        end if
        if (axis == FACES) then
            width = width + 1
        end if

        if (SIZE(coords) /= width) then
           call Driver_abortFlash("Get Coords: coordinates array is wrong size")
        end if
#endif

        if     ((edge == LEFT_EDGE) .OR. (edge == FACES))  then
            shift = 1.0
        else if (edge == CENTER) then
            shift = 0.5
        else if (edge == RIGHT_EDGE) then
            shift = 0.0
#ifdef DEBUG_GRID
        else
            call Driver_abortFlash('[Grid_getCellCoods] invalid edge')
#endif
        end if

        if (guardcell) then
            shift = shift - NGUARD
        end if

        associate (x0   => amrex_problo(axis), &
                   x_lo => this%limits(LOW, axis) - 1, &
                   dx   => amrex_geom(this%level - 1)%dx(axis))
            do i = 1, SIZE(coords)
                coords(i) = x0 + (x_lo + i - shift) * dx
            end do
        end associate
    end subroutine coordinates

    function enclosingBlock(this)
        use amrex_fort_module,    ONLY : wp => amrex_real
        use amrex_box_module,     ONLY : amrex_box

        use gr_physicalMultifabs, ONLY : unk

        class(flash_tile_t), intent(IN)  :: this
        type(flash_tile_t)               :: enclosingBlock

        type(amrex_box) :: parent_box
        integer         :: idx(1:MDIM+1)

        real(wp), pointer :: parent_dataPtr(:, :, :, :)

        enclosingBlock%level      = this%level
        enclosingBlock%grid_index = this%grid_index
        enclosingBlock%tile_index = 0

        associate(level => this%level - 1, &
                  gId   => this%grid_index)
            ! The boxes in the boxarray describe the interiors of
            ! our blocks
            parent_box = unk(level)%ba%get_box(gId)
            enclosingBlock%limits(:, :) = 1
            enclosingBlock%limits(LOW,  1:NDIM) = parent_box%lo(1:NDIM) + 1
            enclosingBlock%limits(HIGH, 1:NDIM) = parent_box%hi(1:NDIM) + 1

            call parent_box%grow(NGUARD)
            enclosingBlock%limitsGC(:, :) = 1
            enclosingBlock%limitsGC(LOW,  1:NDIM) = parent_box%lo(1:NDIM) + 1
            enclosingBlock%limitsGC(HIGH, 1:NDIM) = parent_box%hi(1:NDIM) + 1

            enclosingBlock%blkLimitsGC(:, :) = 1
            enclosingBlock%blkLimitsGC(LOW,  1:NDIM) = parent_box%lo(1:NDIM) + 1
            enclosingBlock%blkLimitsGC(HIGH, 1:NDIM) = parent_box%hi(1:NDIM) + 1
        end associate
    end function enclosingBlock

    ! DEV: If the client code requests a pointer to data that is not 
    ! included in the problem, this routine will return a null pointer
    ! without indicating an error.
    !
    ! This gives the client code the possibility to use either preprocessor
    ! checks to avoid calling this routine needlessly or to do runtime checks
    ! of pointers.
    subroutine getDataPtr(this, dataPtr, gridDataStruct)
        use amrex_fort_module,      ONLY : wp => amrex_real

        use gr_physicalMultifabs,   ONLY : unk, &
                                           gr_scratchCtr, &
                                           facevarx, facevary, facevarz, &
                                           fluxes

        class(flash_tile_t), intent(IN),  target   :: this
        real(wp),                         pointer  :: dataPtr(:, :, :, :)
        integer,             intent(IN)            :: gridDataStruct

        integer :: lo(1:MDIM)

        ! Avoid possible memory leaks
        if (associated(dataPtr)) then
            call Driver_abortFlash("[getDataPtr] Given data pointer must be NULL")
        end if

        lo = this%blkLimitsGC(LOW, :)
        ! These multifabs are hardwired at creation so that the FAB data only
        ! exists for the block interiors
        if (     (gridDataStruct == SCRATCH_CTR) .OR. (gridDataStruct == FLUXX) & 
            .OR. (gridDataStruct == FLUXY)       .OR. (gridDataStruct == FLUXZ)) then
           lo(1:NDIM) = lo(1:NDIM) + NGUARD
        end if
  
        ! Multifab arrays use 0-based level index set (AMReX) instead of 
        ! 1-based set (FLASH/block)
        associate (ilev => this%level - 1, &
                   igrd => this%grid_index)
          select case (gridDataStruct)
          case(CENTER)
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => unk     (ilev)%dataptr(igrd)
          case(FACEX)
#if NFACE_VARS > 0
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevarx(ilev)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(FACEY)
#if NFACE_VARS > 0 && NDIM >= 2
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevary(ilev)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(FACEZ)
#if NFACE_VARS > 0 && NDIM == 3
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevarz(ilev)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(FLUXX)
#if NFLUXES > 0
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, IAXIS)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(FLUXY)
#if NFLUXES > 0 && NDIM >= 2
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, JAXIS)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(FLUXZ)
#if NFLUXES > 0 && NDIM == 3
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, KAXIS)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(SCRATCH_CTR)
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => gr_scratchCtr(ilev)%dataptr(igrd)
          case DEFAULT
              call Driver_abortFlash("[getDataPtr] Unknown grid data structure")
          end select
        end associate
    end subroutine getDataPtr

    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        use amrex_fort_module, ONLY : wp => amrex_real

        class(flash_tile_t), intent(IN)            :: this
        real(wp),            intent(OUT), pointer  :: dataPtr(:, :, :, :)
        integer,             intent(IN)            :: gridDataStruct

        nullify(dataPtr)
    end subroutine releaseDataPtr

end module flash_tile

