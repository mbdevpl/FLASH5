!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getBlkPtr
!!
!! NAME
!!  Grid_getBlkPtr
!!
!! SYNOPSIS
!!
!!  Grid_getBlkPtr(block_metadata_t(IN)   :: block,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!!  Gets a pointer to a single block of simulation data from the
!!  specified Grid data structure. The block includes guard cells.
!!  If the optional argument "gridDataStructure" is not specified,
!!  it returns a block from cell centered data structure.
!!
!!  When using Paramesh 4 in NO_PERMANENT_GUARDCELLS mode, it is important to
!!  release the block pointer for a block before getting it for another block.
!!  For example if pointer to block 1 is not yet released and the user
!!  tries to get a pointer to block 2, the routine will abort.
!!
!! ARGUMENTS 
!!
!!  block : derived type containing metadata for block whose data we need to
!!          access
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space facex variables
!!                   SCRATCH_FACEY scratch space facey variables
!!                   SCRATCH_FACEZ scratch space facez variables
!!
!!
!!
!! NOTES
!!
!!  Grid_getBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_releaseBlkPtr when you are finished with it!
!!
!!***

!!REORDER(4): dataPtr

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getBlkPtr(blockID, dataPtr, gridDataStruct)
  use Driver_interface, ONLY : Driver_abortFlash
  
  implicit none

  integer, intent(in) :: blockID
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, optional,intent(in) :: gridDataStruct
 
  call Driver_abortFlash("AMReX does *not* deal in block IDs")
end subroutine Grid_getBlkPtr

subroutine Grid_getBlkPtr_Itor(block, dataPtr, gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use block_metadata, ONLY : block_metadata_t
  ! DEVNOTE:  Once we have mesh initialization done, we will need to get access
  ! to AMReX owned physical data here
!  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
!  use gr_specificData, ONLY : scratch,scratch_ctr,&
!       scratch_facevarx,scratch_facevary,scratch_facevarz

  implicit none

  type(block_metadata_t), intent(in)            :: block
  real,                   intent(out), pointer  :: dataPtr(:, :, :, :)
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

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

  associate (lo   => block%limitsGC(LOW, :), &
             ilev => block%level, &
             igrd => block%grid_index)
!#ifdef INDEXREORDER
    ! TODO: Code up reordered indices once the normal way is finalized
!#else
    ! TODO: This is just an idea of how to get data.  Get data from the 
    ! AMReX Multifabs correctly once the mesh initialization is done.
    select case (gds)
    case(CENTER)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => unk(ilev)%dataptr(igrd)
    case(FACEX)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx(ilev)%dataptr(igrd)
    case(FACEY)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary(ilev)%dataptr(igrd)
    case(FACEZ)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz(ilev)%dataptr(igrd)
    case(SCRATCH)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch(ilev)%dataptr(igrd)
    case(SCRATCH_CTR)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_ctr(ilev)%dataptr(igrd)
    case(SCRATCH_FACEX)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarx(ilev)%dataptr(igrd)
    case(SCRATCH_FACEY)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevary(ilev)%dataptr(igrd)
    case(SCRATCH_FACEZ)
       dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarz(ilev)%dataptr(igrd)
    case(WORK)
       call Driver_abortFlash("work array cannot be got as pointer")
    case DEFAULT
       print *, 'TRIED TO GET SOMETHING OTHER THAN UNK OR SCRATCH OR FACE[XYZ]. NOT YET.'
    end select
!#endif
  end associate
end subroutine Grid_getBlkPtr_Itor

