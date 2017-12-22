!!****if* source/Grid/GridMain/Grid_getBlkPtr
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
!!FOR FUTURE: Add REORDER for unk, facevar[xyz]1, etc.?

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getBlkPtr_desc(block, dataPtr, gridDataStruct,localFlag)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use gr_specificData, ONLY : scratch,scratch_ctr,&
       scratch_facevarx,scratch_facevary,scratch_facevarz
  use block_metadata, ONLY : block_metadata_t
  use workspace, ONLY : work
#ifdef FL_NON_PERMANENT_GUARDCELLS
  use physicaldata, ONLY: unk1, facevarx1, facevary1, facevarz1,&
                          gcell_on_cc, gcell_on_fc
  use paramesh_interfaces, ONLY : amr_1blk_guardcell
  use Grid_data, ONLY : gr_meshMe, gr_blkPtrRefCount, gr_lastBlkPtrGotten, &
       gr_blkPtrRefCount_fc, gr_lastBlkPtrGotten_fc,gr_ccMask,gr_fcMask
#endif 
  use amrex_multifab_module
  use gr_amrextData, ONLY : gr_amrextUnkMFs

  implicit none
  type(block_metadata_t), intent(in),TARGET :: block
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, optional,intent(in) :: gridDataStruct
  logical,optional, intent(in) :: localFlag
  
  type(amrex_multifab), POINTER :: mf
  integer :: gds, blkPtrRefCount, lastBlkPtrGotten
  logical :: validGridDataStruct
  integer,pointer,dimension(:) :: loUse
  integer :: blockID

#ifdef FL_NON_PERMANENT_GUARDCELLS
  integer :: idest, iopt, nlayers, icoord
  logical :: lcc, lfc, lec, lnc, l_srl_only, ldiag
  logical,dimension(NUNK_VARS) :: save_ccMask
#if NFACE_VARS >0
  logical, dimension(3,NFACE_VARS) :: save_fcMask
#endif
#endif

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

#ifdef FL_NON_PERMANENT_GUARDCELLS
  if (gds .eq. CENTER .or. gds .eq. FACEX .or. gds .eq. FACEY .or. gds .eq. FACEZ) then
     idest = 1
     iopt = 1
     nlayers = NGUARD
     if (gds .eq. FACEX .or. gds .eq. FACEY .or. gds .eq. FACEZ) then
        blkPtrRefCount = gr_blkPtrRefCount_fc
        lastBlkPtrGotten = gr_lastBlkPtrGotten_fc
#if NFACE_VARS>0
        save_fcMask=gcell_on_fc
        gcell_on_fc=gr_fcMask
#endif
        lcc = .false.
        lfc = .true.
     else
        blkPtrRefCount = gr_blkPtrRefCount
        lastBlkPtrGotten = gr_lastBlkPtrGotten
        save_ccMask=gcell_on_cc
        gcell_on_cc=gr_ccMask
        lcc = .true.
        lfc = .false.
     end if
     lec = .false.
     lnc = .false.
     l_srl_only = .false.
     icoord = 0
     ldiag = .true.

     if (blkPtrRefCount .ne. 0 ) then 
        if (blockId .ne. lastBlkPtrGotten) then
           call Driver_abortFlash("Grid_getBlkPtr: you can't get another pointer while one's in use, " // &
                "unless it's to the block that's in use")
        end if
     else
        blkPtrRefCount = 0
        lastBlkPtrGotten = blockId

        call amr_1blk_guardcell(gr_meshMe,iopt,nlayers,blockId,gr_meshMe, &
             lcc,lfc,lec,lnc, &
             l_srl_only,icoord,ldiag)
     end if

     blkPtrRefCount = blkPtrRefCount + 1
     if(gds==CENTER) then
        gcell_on_cc=save_ccMask
     else
#if NFACE_VARS>0
        gcell_on_fc=save_fcMask
#endif
     end if
  end if
#endif 
  !  end of #ifdef FL_NON_PERMANENT_GUARDCELLS
     
  loUse => block%limitsGC(LOW, :)
  if (present(localFlag)) then
     if (localFlag) loUse => block%localLimitsGC(LOW, :)
  end if

     associate (lo => loUse)
#ifdef INDEXREORDER
        select case (gds)
#ifndef FL_NON_PERMANENT_GUARDCELLS
        case(CENTER)
           if (block%grid_index < 0) then
!!$              print*,'Grid_getBlkPtr_desc: no valid grid_index!!!!!' ! An annoying warning that we are called with a synthetic block descriptor...
              dataPtr(1:, lo(1):, lo(2):, lo(3):) => block%fp
           else
              mf => gr_amrextUnkMFs(block%level-1)
!!$              print*,'Grid_getBlkPtr: (level,grid_index)=',block%level,block%grid_index
              dataPtr(1:, lo(1):, lo(2):, lo(3):) => mf%dataPtr(block%grid_index)
           end if

!!$        case(FACEX)
!!$           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx(ilev)%dataptr(igrd)
!!$        case(FACEY)
!!$           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary(ilev)%dataptr(igrd)
!!$        case(FACEZ)
!!$           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz(ilev)%dataptr(igrd)
#else
        !  #ifndef FL_NON_PERMANENT_GUARDCELLS ...
        case(CENTER)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => unk1(:,:,:,:,idest)
        case(FACEX)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx1(:,:,:,:,idest)
        case(FACEY)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary1(:,:,:,:,idest)
        case(FACEZ)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz1(:,:,:,:,idest)
#endif
        !  end of #ifdef FL_NON_PERMANENT_GUARDCELLS
        case(SCRATCH)
           dataPtr(lo(1):, lo(2):, lo(3):, 1:) => scratch(:,:,:,:,blockid)
        case(SCRATCH_CTR)
           dataPtr(lo(1):, lo(2):, lo(3):, 1:) => scratch_ctr(:,:,:,:,blockid)           
        case(SCRATCH_FACEX)
           dataPtr(lo(1):, lo(2):, lo(3):, 1:) => scratch_facevarx(:,:,:,:,blockid)           
        case(SCRATCH_FACEY)
           dataPtr(lo(1):, lo(2):, lo(3):, 1:) => scratch_facevary(:,:,:,:,blockid)           
        case(SCRATCH_FACEZ)
           dataPtr(lo(1):, lo(2):, lo(3):, 1:) => scratch_facevarz(:,:,:,:,blockid)           
        case DEFAULT
           print *, 'TRIED TO GET SOMETHING OTHER THAN UNK OR SCRATCH OR FACE[XYZ]. NOT YET.'
        case(WORK)
           call Driver_abortFlash("work array cannot be got as pointer")
        end select
#else
        select case (gds)
#ifndef FL_NON_PERMANENT_GUARDCELLS
        case(CENTER)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => block%fp
!!$        case(FACEX)
!!$           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx(ilev)%dataptr(igrd)
!!$        case(FACEY)
!!$           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary(ilev)%dataptr(igrd)
!!$        case(FACEZ)
!!$           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz(ilev)%dataptr(igrd)
#else
        !  #ifndef FL_NON_PERMANENT_GUARDCELLS ...
        case(CENTER)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => unk1(:,:,:,:,idest)
        case(FACEX)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx1(:,:,:,:,idest)
        case(FACEY)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary1(:,:,:,:,idest)
        case(FACEZ)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz1(:,:,:,:,idest)
#endif
        !  end of #ifdef FL_NON_PERMANENT_GUARDCELLS
        case(SCRATCH)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch(:,:,:,:,blockid)
        case(SCRATCH_CTR)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_ctr(:,:,:,:,blockid)           
        case(SCRATCH_FACEX)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarx(:,:,:,:,blockid)           
        case(SCRATCH_FACEY)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevary(:,:,:,:,blockid)           
        case(SCRATCH_FACEZ)
           dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarz(:,:,:,:,blockid)           
        case DEFAULT
           print *, 'TRIED TO GET SOMETHING OTHER THAN UNK OR SCRATCH OR FACE[XYZ]. NOT YET.'
        case(WORK)
           call Driver_abortFlash("work array cannot be got as pointer")
        end select
#endif
     end associate

#ifdef FL_NON_PERMANENT_GUARDCELLS
  if (gds .eq. FACEX .or. gds .eq. FACEY .or. gds .eq. FACEZ) then
     gr_blkPtrRefCount_fc = blkPtrRefCount
     gr_lastBlkPtrGotten_fc = lastBlkPtrGotten
  else if (gds .eq. CENTER) then
     gr_blkPtrRefCount = blkPtrRefCount
     gr_lastBlkPtrGotten = lastBlkPtrGotten
  end if
#endif

  return
end subroutine Grid_getBlkPtr_desc

