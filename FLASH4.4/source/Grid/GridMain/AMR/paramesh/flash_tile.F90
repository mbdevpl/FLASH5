!!****ih* source/Grid/GridMain/AMR/paramesh/flash_tile
!!
!!
!!
!!****

!!REORDER(4): dataPtr

#include "constants.h"
#include "Flash.h"

module flash_tile
    implicit none

    private

    type, public :: flash_tile_t 
        integer :: id
        integer :: cid(MDIM)
        integer :: stride(MDIM)
        integer :: level
        integer :: limits(LOW:HIGH, MDIM)
        integer :: grownLimits(LOW:HIGH, MDIM)
        integer :: blkLimitsGC(LOW:HIGH, MDIM)
    contains
        procedure, public :: deltas
        procedure, public :: coordinates
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
    end type flash_tile_t

contains

    subroutine deltas(this, dx)
        use Grid_data, ONLY: gr_delta

        class(flash_tile_t), intent(IN)  :: this
        real,                intent(OUT) :: dx(1:MDIM)

        dx(1:MDIM) = gr_delta(1:MDIM, this%level)
    end subroutine deltas

    subroutine coordinates(this, axis, edge, region, coords)
        use Grid_data,        ONLY : gr_globalDomain, &
                                     gr_delta
        use Driver_interface, ONLY : Driver_abortFlash

        class(flash_tile_t), intent(IN)  :: this
        integer,             intent(IN)  :: axis
        integer,             intent(IN)  :: edge
        integer,             intent(IN)  :: region 
        real,                intent(OUT) :: coords(:)

        real    :: shift
        integer :: nCells
        integer :: lo
        integer :: i

#ifdef DEBUG_GRID
        print*,' get coordinates', axis, this%id, edge, guardcell,size
        if ((this%id < 1) .OR. (this%id > MAXBLOCKS)) then
           call Driver_abortFlash("Grid_getCellCoords :invalid blockID ")
        end if
        if(.not.((edge==LEFT_EDGE).or.(edge==RIGHT_EDGE).or.(edge==FACES).or.&
             &(edge==CENTER))) then
           call Driver_abortFlash("Get Coords : invalid edge specification, must be LEFT_EDGE &
                RIGHT_EDGE, or CENTER")
        end if

        if(.not.((axis==IAXIS).or.(axis==JAXIS).or.(axis==KAXIS))) then
           call Driver_abortFlash("Get Coords : invalid axis, must be IAXIS, JAXIS or KAXIS ")
        end if
#endif

        ! No tiling with Paramesh so that GROWN_TILE and TILE_AND_HALO
        ! yield the same result - the block and its GC halo
        if       (region == TILE) then
           lo = this%limits(LOW, axis)
           nCells = this%limits(HIGH, axis) - this%limits(LOW, axis) + 1
        else if ((region == GROWN_TILE) .OR. (region == TILE_AND_HALO)) then
           lo = this%blkLimitsGC(LOW, axis)
           nCells = this%blkLimitsGC(HIGH, axis) - this%blkLimitsGC(LOW, axis) + 1
        else
          coords(:) = 0.0
          call Driver_abortFlash("[coordinates] Unknown region")
        end if

        if      (edge == FACES)  then
            shift = 2.0
            nCells = nCells + 1
        else if (edge == LEFT_EDGE) then
            shift = 2.0
        else if (edge == CENTER) then
            shift = 1.5
        else if (edge == RIGHT_EDGE) then
            shift = 1.0
        else
            call Driver_abortFlash('[coordinates] Invalid edge')
        end if

        associate (x0 => gr_globalDomain(LOW, axis), &
                   dx => gr_delta(axis, this%level))
            do i = 1, nCells
                coords(i) = x0 + (lo + i - shift) * dx
            end do
        end associate
    end subroutine coordinates

    subroutine getDataPtr(this, dataPtr, gridDataStruct)
        use physicaldata,    ONLY : unk, &
                                    facevarx, facevary, facevarz
        use gr_specificData, ONLY : scratch, &
                                    scratch_ctr, &
                                    scratch_facevarx, &
                                    scratch_facevary, & 
                                    scratch_facevarz, &
                                    gr_flxx, gr_flxy, gr_flxz

       class(flash_tile_t), intent(IN),  target   :: this
       real,                             pointer  :: dataPtr(:, :, :, :)
       integer,             intent(IN)            :: gridDataStruct

#ifdef DEBUG_GRID
       validGridDataStruct = .false.
       validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_CTR).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_FACEX).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_FACEY).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_FACEZ).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
       if(.NOT. validGridDataStruct) then
          print *, "Grid_getBlkPtr: gridDataStruct set to improper value"
          print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ," // &
               "WORK or SCRATCH (defined in constants.h)"
          call Driver_abortFlash("gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,SCRATCH (see constants.h)")
       end if

       if ((this%id < 1) .OR. (this%id > MAXBLOCKS)) then
          print *, 'Grid_getBlkPtr:  invalid blockid ',this%id
          call Driver_abortFlash("[getDataPtr] invalid blockid ")
       end if
#endif

       ! Avoid possible memory leaks
       if (associated(dataPtr)) then
           call Driver_abortFlash("[getDataPtr] Given data pointer must be NULL")
       end if

       ! Note that loFl is presently set under the assumption that tiling
       ! is *not* implemented with Paramesh
       associate(blockID => this%id, &
                 lo      => this%blkLimitsGC(LOW, :), &
                 loFl    => this%limits(LOW, :))
          select case (gridDataStruct)
          case(CENTER)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => unk(:,:,:,:,blockID)
          case(FACEX)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx(:,:,:,:,blockID)
          case(FACEY)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary(:,:,:,:,blockID)
          case(FACEZ)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz(:,:,:,:,blockID)
          case(FLUXX)
             dataPtr(1:, loFl(1):, loFl(2):, loFl(3):) => gr_flxx(:,:,:,:,blockID)
          case(FLUXY)
             dataPtr(1:, loFl(1):, loFl(2):, loFl(3):) => gr_flxy(:,:,:,:,blockID)
          case(FLUXZ)
             dataPtr(1:, loFl(1):, loFl(2):, loFl(3):) => gr_flxz(:,:,:,:,blockID)
          case(SCRATCH)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch(:,:,:,:,blockID)
          case(SCRATCH_CTR)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_ctr(:,:,:,:,blockID)           
          case(SCRATCH_FACEX)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarx(:,:,:,:,blockID)           
          case(SCRATCH_FACEY)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevary(:,:,:,:,blockID)           
          case(SCRATCH_FACEZ)
             dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarz(:,:,:,:,blockID)           
          case DEFAULT
             print *, 'TRIED TO GET SOMETHING OTHER THAN UNK OR SCRATCH OR FACE[XYZ]. NOT YET.'
          case(WORK)
             call Driver_abortFlash("work array cannot be got as pointer")
          end select
       end associate
    end subroutine getDataPtr

    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        class(flash_tile_t), intent(IN)            :: this
        real,                             pointer  :: dataPtr(:, :, :, :)
        integer,             intent(IN)            :: gridDataStruct

        nullify(dataPtr)
    end subroutine releaseDataPtr

end module flash_tile

