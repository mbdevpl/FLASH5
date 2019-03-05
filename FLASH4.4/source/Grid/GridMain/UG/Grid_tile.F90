!!****ih* source/Grid/GridMain/UG/Grid_tile
!!
!!
!!****

!!REORDER(5): unk, facevar[xyz]
!!REORDER(4): dataPtr

#include "constants.h"
#include "Flash.h"

module Grid_tile
    implicit none

    private

    type, public :: Grid_tile_t
        integer, public :: id
        integer, public :: cid(MDIM)
        integer, public :: stride(MDIM)
        integer, public :: level
        integer, public :: limits(LOW:HIGH, MDIM)
        integer, public :: grownLimits(LOW:HIGH, MDIM)
        integer, public :: blkLimitsGC(LOW:HIGH, MDIM)
    contains
        procedure, public :: deltas
        procedure, public :: boundBox
        procedure, public :: physicalSize
        procedure, public :: faceBCs
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
        procedure, public :: enclosingBlock
    end type Grid_tile_t

contains

    subroutine deltas(this, dx)
        use Grid_data, ONLY : gr_delta

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: dx(1:MDIM)

        dx = gr_delta(:, 1)
    end subroutine deltas

    subroutine boundBox(this, box)
        use Grid_data, ONLY : gr_iCoords, gr_jCoords, gr_kCoords
        use Grid_data, ONLY : gr_ilo, gr_ihi, &
                              gr_jlo, gr_jhi, &
                              gr_klo, gr_khi

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: box(LOW:HIGH, 1:MDIM)

        box(:, :) = 0.0
   
        box(LOW,  IAXIS) = gr_iCoords(LEFT_EDGE,  gr_ilo, 1)
        box(HIGH, IAXIS) = gr_iCoords(RIGHT_EDGE, gr_ihi, 1)
#if NDIM > 1
        box(LOW,  JAXIS) = gr_jCoords(LEFT_EDGE,  gr_jlo, 1)
        box(HIGH, JAXIS) = gr_jCoords(RIGHT_EDGE, gr_jhi, 1)
#endif
#if NDIM > 2
        box(LOW,  KAXIS) = gr_kCoords(LEFT_EDGE,  gr_klo, 1)
        box(HIGH, KAXIS) = gr_kCoords(RIGHT_EDGE, gr_khi, 1)
#endif
    end subroutine boundBox 

    subroutine physicalSize(this, tileSize) 
        use Grid_data, ONLY : gr_iCoords, gr_jCoords, gr_kCoords
        use Grid_data, ONLY : gr_ilo, gr_ihi, &
                              gr_jlo, gr_jhi, &
                              gr_klo, gr_khi

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: tileSize(1:MDIM) 

        tileSize(:) = 0.0
        tileSize(IAXIS) =   gr_iCoords(RIGHT_EDGE, gr_ihi, 1) &
                          - gr_iCoords(LEFT_EDGE,  gr_ilo, 1)
#if NDIM > 1
        tileSize(JAXIS) =   gr_jCoords(RIGHT_EDGE, gr_jhi, 1) &
                          - gr_jCoords(LEFT_EDGE,  gr_jlo, 1)
#endif
#if NDIM > 2
        tileSize(KAXIS) =   gr_kCoords(RIGHT_EDGE, gr_khi, 1) &
                          - gr_kCoords(LEFT_EDGE,  gr_klo, 1)
#endif
    end subroutine physicalSize
    
    subroutine faceBCs(this, faces, onBoundary)
        use Grid_data, ONLY : gr_blkBC

        class(Grid_tile_t), intent(IN)            :: this
        integer,            intent(OUT)           :: faces(LOW:HIGH, 1:MDIM)
        integer,            intent(OUT), optional :: onBoundary(LOW:HIGH, 1:MDIM)

        faces = gr_blkBC
        where (faces == PERIODIC)
           faces = NOT_BOUNDARY
        end where 

        if (present(onBoundary)) then     
           onBoundary = gr_blkBC
        end if
    end subroutine faceBCs

    function enclosingBlock(this)
        use Driver_interface, ONLY : Driver_abortFlash

        class(Grid_tile_t), intent(IN) :: this
        type(Grid_tile_t)              :: enclosingBlock
            
        call Driver_abortFlash("[enclosingBlock] not implemented yet")
    end function enclosingBlock

    ! DEV: If the client code requests a pointer to data that is not 
    ! included in the problem, this routine will return a null pointer
    ! without indicating an error.
    !
    ! This gives the client code the possibility to use either preprocessor
    ! checks to avoid calling this routine needlessly or to do runtime checks
    ! of pointers.
    !
    ! DEV: For now, the localFlag parameter might be useful as we pull in more
    !      units from FLASH4.4.  Try to get rid of it along the way.
    subroutine getDataPtr(this, dataPtr, gridDataStruct, localFlag)
        use physicaldata, ONLY : unk, &
                                 facevarx, facevary, facevarz
        use Grid_data,    ONLY : gr_flxx, gr_flxy, gr_flxz, &
                                 gr_ilo,   gr_jlo,   gr_klo, &
                                 gr_iloGc, gr_jloGc, gr_kloGc

        class(Grid_tile_t), intent(IN),  target   :: this
        real,                            pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct
        logical,            intent(IN),  optional :: localFlag

        ! Avoid possible memory leaks
        if (associated(dataPtr)) then
            call Driver_abortFlash("[getDataPtr] Given data pointer must be NULL")
        end if

        if (present(localFlag)) then
            call Driver_abortFlash("[getDataPtr] localFlag not implemented")
        end if

        if (this%level /= 1) then
            call Driver_abortFlash("[getDataPtr] Level must be one")
        end if

        select case (gridDataStruct)
        case(CENTER)
           dataPtr(1:, gr_iloGc:, gr_jloGc:, gr_kloGc:) => unk(:,:,:,:,1)
        case(FACEX)
           dataPtr => facevarx(1:,gr_iloGc:,gr_jloGc:,gr_kloGc:,1)
        case(FACEY)
           dataPtr => facevary(1:,gr_iloGc:,gr_jloGc:,gr_kloGc:,1)
        case(FACEZ)
           dataPtr => facevarz(1:,gr_iloGc:,gr_jloGc:,gr_kloGc:,1)
        case(FLUXX)
           dataPtr(1:, gr_ilo:, gr_jlo:, gr_klo:) => gr_flxx(:,:,:,:)
        case(FLUXY)
           dataPtr(1:, gr_ilo:, gr_jlo:, gr_klo:) => gr_flxy(:,:,:,:)
        case(FLUXZ)
           dataPtr(1:, gr_ilo:, gr_jlo:, gr_klo:) => gr_flxz(:,:,:,:) 
        case DEFAULT
            call Driver_abortFlash("[getDataPtr] Unknown grid data structure")
        end select
    end subroutine getDataPtr

    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        class(Grid_tile_t), intent(IN)            :: this
        real,               intent(OUT), pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct

        nullify(dataPtr)
    end subroutine releaseDataPtr

end module Grid_tile

