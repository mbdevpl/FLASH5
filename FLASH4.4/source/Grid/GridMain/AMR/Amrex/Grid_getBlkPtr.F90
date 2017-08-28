subroutine Grid_getBlkPtr(blockID, blkPtr, gridDataStruct)
  implicit none

  integer, intent(IN)            :: blockID
  real,    intent(OUT), pointer  :: blkPtr(:, :, :, :)
  integer, intent(IN),  optional :: gridDataStruct

  write(*,*) "AMReX does *not* deal in block IDs"
  stop
end subroutine Grid_getBlkPtr

subroutine Grid_getBlkPtr_Itor(block, blkPtr, gridDataStruct)

!#include "constants.h"
!#include "Flash.h"

  use amrex_fort_module, ONLY : wp => amrex_real

!  use Driver_interface, ONLY : Driver_abortFlash
  use physicaldata,      ONLY : unk
  use block_metadata,    ONLY : block_metadata_t

  implicit none

#include "constants.h"

  ! DEV: How to match data types for blkPtr with FLASH?
  type(block_metadata_t), intent(in)            :: block
  real(wp),               intent(out), pointer  :: blkPtr(:, :, :, :)
  integer,                intent(in),  optional :: gridDataStruct

  integer :: gds
  logical :: validGridDataStruct

#ifdef DEBUG_GRID
  if(present(gridDataStruct)) then
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
     
     if(.not.validGridDataStruct) then
        print *, "Grid_getBlkPtr: gridDataStruct set to improper value"
        print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ," // &
             "WORK or SCRATCH (defined in constants.h)"
        call Driver_abortFlash("gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,SCRATCH (see constants.h)")
     end if
  end if
  ! TODO: Convert this into error checking of AMReX metadata
!  if((blockid<1).or.(blockid>MAXBLOCKS)) then
!     print *, 'Grid_getBlkPtr:  invalid blockid ',blockid
!     call Driver_abortFlash("[Grid_getBlkPtr] invalid blockid ")
!  end if
#endif

!  if(present(gridDataStruct)) then
!     gds = gridDataStruct
!  else
!     gds = CENTER
!  end if

  associate (lo   => block%limitsGC(LOW, :), &
             ilev => block%level, &
             igrd => block%grid_index)
!    select case (gds)
!    case(CENTER)
       blkPtr => unk(ilev)%dataptr(igrd)
!    case(FACEX)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx(ilev)%dataptr(igrd)
!    case(FACEY)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => facevary(ilev)%dataptr(igrd)
!    case(FACEZ)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz(ilev)%dataptr(igrd)
!    case(SCRATCH)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => scratch(ilev)%dataptr(igrd)
!    case(SCRATCH_CTR)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_ctr(ilev)%dataptr(igrd)
!    case(SCRATCH_FACEX)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarx(ilev)%dataptr(igrd)
!    case(SCRATCH_FACEY)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevary(ilev)%dataptr(igrd)
!    case(SCRATCH_FACEZ)
!       blkPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarz(ilev)%dataptr(igrd)
!    case(WORK)
!       call Driver_abortFlash("work array cannot be got as pointer")
!    case DEFAULT
!       print *, 'TRIED TO GET SOMETHING OTHER THAN UNK OR SCRATCH OR FACE[XYZ]. NOT YET.'
!    end select
  end associate
end subroutine Grid_getBlkPtr_Itor

