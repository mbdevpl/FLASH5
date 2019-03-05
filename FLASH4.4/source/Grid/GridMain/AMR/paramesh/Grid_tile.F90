!!****ih* source/Grid/GridMain/AMR/paramesh/Grid_tile
!!
!!
!!
!!****

!!REORDER(4): dataPtr

#include "constants.h"
#include "Flash.h"

module Grid_tile
    implicit none

    private

    type, public :: Grid_tile_t 
        integer :: id
        integer :: cid(MDIM)
        integer :: stride(MDIM)
        integer :: level
        integer :: limits(LOW:HIGH, MDIM)
        integer :: grownLimits(LOW:HIGH, MDIM)
        integer :: blkLimitsGC(LOW:HIGH, MDIM)
    contains
        procedure, public :: deltas
        procedure, public :: boundBox
        procedure, public :: physicalSize
        procedure, public :: faceBCs
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
    end type Grid_tile_t

contains

    subroutine deltas(this, dx)
        use Grid_data, ONLY: gr_delta

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: dx(1:MDIM)

        dx(1:MDIM) = gr_delta(1:MDIM, this%level)
    end subroutine deltas

    subroutine boundBox(this, box)
        use tree,             ONLY : bnd_box
        use Driver_interface, ONLY : Driver_abortFlash

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: box(LOW:HIGH, 1:MDIM)
  
        if (this%id <= 0) then
           print *, "blockId = ", this%id
           call Driver_abortFlash("[boundBox] blockId out of bounds")
        end if
  
        box = bnd_box(:, :, this%id)
    end subroutine boundBox 

    subroutine physicalSize(this, tileSize) 
        use tree, ONLY : bsize

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: tileSize(1:MDIM) 
      
        tileSize = bsize(:, this%id)
    end subroutine physicalSize

    subroutine faceBCs(this, faces, onBoundary)
        use tree,         ONLY : bnd_box
        use Grid_data,    ONLY : gr_domainBC, &
                                 gr_globalDomain, &
                                 gr_delta

        class(Grid_tile_t), intent(IN)            :: this
        integer,            intent(OUT)           :: faces(LOW:HIGH, 1:MDIM)
        integer,            intent(OUT), optional :: onBoundary(LOW:HIGH, 1:MDIM)

        real    :: deltas(1:MDIM)
        integer :: axis, face

        deltas(1:MDIM) = gr_delta(1:MDIM, this%level)

        do    axis = 1, MDIM
           do face = LOW, HIGH
              faces(face, axis) = NOT_BOUNDARY
              if (present (onBoundary)) then
                 onBoundary(face,axis) = NOT_BOUNDARY
              end if

              if (almostEqual(bnd_box(face, axis, this%id), &
                              gr_globalDomain(face, axis), &
                              deltas(axis))) then
                 if (gr_domainBC(face, axis) .NE. PERIODIC) &
                      faces(face,axis) = gr_domainBC(face,axis)
                 if (present (onBoundary)) then
                    onBoundary(face,axis) = gr_domainBC(face,axis)
                 end if
              end if

           end do
        end do

    contains
       logical function almostEqual(x, y, dx)
          real, intent(IN) :: x
          real, intent(IN) :: y
          real, intent(IN) :: dx

          almostEqual = (ABS(x-y) <= (0.01 * dx))
       end function almostEqual
    end subroutine faceBCs

    subroutine getDataPtr(this, dataPtr, gridDataStruct, localFlag)
        use physicaldata,    ONLY : unk, &
                                    facevarx, facevary, facevarz
        use gr_specificData, ONLY : scratch, &
                                    scratch_ctr, &
                                    scratch_facevarx, &
                                    scratch_facevary, & 
                                    scratch_facevarz, &
                                    gr_flxx, gr_flxy, gr_flxz

       class(Grid_tile_t), intent(IN), target   :: this
       real,                           pointer  :: dataPtr(:, :, :, :)
       integer,            intent(IN)           :: gridDataStruct
       logical,            intent(IN), optional :: localFlag

       integer :: lo(1:MDIM)

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

       lo = this%blkLimitsGC(LOW, :)

       ! These multifabs are hardwired at creation so that the FAB data only
       ! exists for the block interiors
       if (     (gridDataStruct == FLUXX) & 
           .OR. (gridDataStruct == FLUXY) &
           .OR. (gridDataStruct == FLUXZ)) then
          lo(1:NDIM) = lo(1:NDIM) + NGUARD
          if (present(localFlag)) then
              if (localFlag) then
                  lo(1:NDIM) = NGUARD + 1
              end if
          end if
       else if (present(localFlag)) then
           if (localFlag) then
               lo(:) = 1
           end if
       end if

       select case (gridDataStruct)
       case(CENTER)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => unk(:,:,:,:,this%id)
       case(FACEX)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx(:,:,:,:,this%id)
       case(FACEY)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary(:,:,:,:,this%id)
       case(FACEZ)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz(:,:,:,:,this%id)
       case(FLUXX)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => gr_flxx(:,:,:,:,this%id)
       case(FLUXY)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => gr_flxy(:,:,:,:,this%id)
       case(FLUXZ)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => gr_flxz(:,:,:,:,this%id)
       case(SCRATCH)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch(:,:,:,:,this%id)
       case(SCRATCH_CTR)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_ctr(:,:,:,:,this%id)           
       case(SCRATCH_FACEX)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarx(:,:,:,:,this%id)           
       case(SCRATCH_FACEY)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevary(:,:,:,:,this%id)           
       case(SCRATCH_FACEZ)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarz(:,:,:,:,this%id)           
       case DEFAULT
          print *, 'TRIED TO GET SOMETHING OTHER THAN UNK OR SCRATCH OR FACE[XYZ]. NOT YET.'
       case(WORK)
          call Driver_abortFlash("work array cannot be got as pointer")
       end select
    end subroutine getDataPtr

    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        class(Grid_tile_t), intent(IN)            :: this
        real,                            pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct

        nullify(dataPtr)
    end subroutine releaseDataPtr

end module Grid_tile

