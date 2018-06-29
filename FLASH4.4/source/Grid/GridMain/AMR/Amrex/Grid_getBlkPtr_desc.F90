!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getBlkPtr
!!
!! NAME
!!  Grid_getBlkPtr
!!
!! SYNOPSIS
!!  Grid_getBlkPtr(block_metadata_t(IN)   :: block,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Gets a pointer to a single block of simulation data from the
!!  specified Grid data structure. The block includes guard cells.
!!  If the optional argument "gridDataStructure" is not specified,
!!  it returns a block from cell centered data structure.
!!
!! ARGUMENTS 
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
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   If not given, then the routine assumes CENTER.
!!
!! NOTES
!!  Grid_getBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_releaseBlkPtr when you are finished with it!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Flash.h"
#include "constants.h"

subroutine Grid_getBlkPtr_desc(block, dataPtr, gridDataStruct,localFlag)
  use amrex_fort_module,      ONLY : wp => amrex_real
  
  use gr_physicalMultifabs,   ONLY : unk, &
                                     gr_scratchCtr, &
                                     facevarx, facevary, facevarz
  use block_metadata,         ONLY : block_metadata_t
  
  use Driver_interface,       ONLY : Driver_abortFlash

  implicit none

  ! DEV: How to match data types for dataPtr with FLASH?
  type(block_metadata_t), intent(in), target    :: block
  real(wp),               intent(out), pointer  :: dataPtr(:, :, :, :)
  integer,                intent(in),  optional :: gridDataStruct
  logical,      optional, intent(in) :: localFlag

  integer :: gds
  logical :: validGridDataStruct
  integer,pointer,dimension(:) :: loUse

#ifdef DEBUG_GRID
  if(present(gridDataStruct)) then
     validGridDataStruct = .false.
     validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_CTR).or.validGridDataStruct

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

  loUse => block%limitsGC(LOW, :)
  if (present(localFlag)) then
     if (localFlag) loUse => block%localLimitsGC(LOW, :)
  end if

  if (gds==SCRATCH_CTR) then    !For now, the assumption that
     !SCRATCH_CTR ist stored without guard cells is hardwired here.
     loUse => block%limits(LOW, :)
     if (present(localFlag)) then
        if (localFlag) loUse => block%localLimits(LOW, :)
     end if
  end if

  ! Multifab arrays use 0-based level index set (AMReX) instead of 
  ! 1-based set (FLASH/block)
  associate (lo   => loUse, &
             ilev => block%level - 1, &
             igrd => block%grid_index)
    select case (gds)
    case(CENTER)
       dataPtr(lo(1):, lo(2):, lo(3):, 1:) => unk     (ilev)%dataptr(igrd)
#if NFACE_VARS > 0
    case(FACEX)
       dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevarx(ilev)%dataptr(igrd)
#if NDIM >= 2
    case(FACEY)
       dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevary(ilev)%dataptr(igrd)
#endif
#if NDIM == 3
    case(FACEZ)
       dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevarz(ilev)%dataptr(igrd)
#endif
#endif
    case(SCRATCH_CTR)
       dataPtr(lo(1):, lo(2):, lo(3):, 1:) => gr_scratchCtr(ilev)%dataptr(igrd)
    case DEFAULT
        call Driver_abortFlash("[Grid_getBlkPtr_desc] Unknown grid data structure")
    end select
  end associate
end subroutine Grid_getBlkPtr_desc

