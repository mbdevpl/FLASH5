!!****if* source/Grid/GridBoundaryConditions/gr_bcApplyToOneFace
!!
!! NAME
!!  gr_bcApplyToOneFace
!!
!! SYNOPSIS
!!
!!  call gr_bcApplyToOneFace(integer(IN)     :: axis,
!!                      integer(IN)          :: bcType,
!!                      integer(IN)          :: gridDataStruct,
!!                      integer(IN)          :: varCount,
!!                      integer(IN)          :: regionType(MDIM)
!!                      block_metadata_t(IN) :: blockDesc,
!!                      integer(IN)          :: idest)
!!  
!! DESCRIPTION 
!!
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    bcType         - the type of boundary condition
!!    gridDataStruct - In PM3 and PM4 it can have values (CENTER,FACEX,FACEY,FACEZ,WORK),
!!                     in UG (CENTER,FACEX,FACEY,FACEZ), and in PM2 (CENTER,WORK).
!!    varCount       - the number of variable in the data structure specified in gridDataStruct
!!    regionType     - The part of the block that is involved in the boundary condition. This integer
!!                     array can have values (LEFT_EDGE, CENTER, RIGHT_EDGE, WHOLE_VEC and NO_VEC)
!!                     for each of the three dimensions of the physical data.
!!                     LEFT_EDGE implies guard cells along lower face. If this value is specified
!!                     for the axis along which the BC are being applies, that means that we 
!!                     are applying BC to the lowerface, if it is one of the other axes, then
!!                     the BC are being applied to one of the corners. Same is true of RIGHT_EDGE
!!                     except that implies the upperface. CENTER, WHOLE_VEC and NO_VEC values
!!                     are valid only for dimensions other than the one specified in "axis". 
!!                     CENTER implies only the interior cells, WHOLE_VEC implies all cells in 
!!                     the block and NO_VEC implies that the correspoding dimension is not
!!                     a part of the region. Normally this value is most likely to be used
!!                     along KAXIS in a 2D problems, and JAXIS and KAXIS in a 1D problem
!!    blockDesc      - Derived type that encapsulates metadata that uniquely
!!                     characterizes local block to be operated on
!!    idest          - this is used in Paramesh 4, where boundary condition handling is applied
!!                     not to a block's solution data in their permanent location (named UNK, etc.),
!!                     but to a buffer that contain temporary copies of a few (normally, two) blocks'
!!                     solution data (named UNK1, etc.). The dummy argument idest selects a buffer
!!                     slot in that case. It is ignored in other Grid implementations.
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_bcApplyToOneFace(axis,bcType,gridDataStruct,varCount,&
                               regionType,blockDesc,idest)

  use Grid_interface, ONLY : Grid_bcApplyToRegion, &
       Grid_bcApplyToRegionSpecialized
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_bcGetRegion, gr_bcPutRegion
  use gr_hgInterface, ONLY : gr_hgMapBcType !!, gr_hg_amr_1blk_bcset_work
  use block_metadata, ONLY : block_metadata_t
  implicit none
 
  integer, intent(in) :: axis,bcType,gridDataStruct,varCount,idest
  integer,dimension(MDIM),intent(IN) :: regionType
  type(block_metadata_t), intent(IN) :: blockDesc
  integer :: face,guard
  integer,dimension(LOW:HIGH,MDIM) :: endPoints
  integer,dimension(REGION_DIM) :: regionSize
  real,allocatable,dimension(:,:,:,:) :: regionData
  integer,dimension(MDIM-1) :: nextDir
  integer :: i,j,k,n,var,sb,se
  integer :: testBcType

  integer,dimension(LOW:HIGH,MDIM) :: blkLimitsGC,blkLimits

  logical :: isFaceVarNormalDir,applied
  logical,allocatable,dimension(:):: mask
!!  real,allocatable,dimension(:) :: cellCenterSweepCoord, secondCoord,thirdCoord

  if(bcType==PERIODIC) return

  if (gridDataStruct==CENTER_FACES .OR. gridDataStruct==FACES) then
     call gr_bcApplyToOneFaceAllGds(axis,bcType,gridDataStruct,varCount,&
          regionType,blockDesc,idest)
     return                  ! DONE - RETURN TO CALLER FROM HERE!
  end if


  if(regionType(axis)==LEFT_EDGE) then
     face=LOW
  else if(regionType(axis)==RIGHT_EDGE) then
     face=HIGH
  else
     call Driver_abortFlash("[gr_bcApplyToOneFace] type along BC dir must be LEFT_EDGE or RIGHT_EDGE")
  end if

#if !defined(FLASH_GRID_PARAMESH2) && !defined(FLASH_GRID_AMREX)
  if(gridDataStruct == WORK) then
     call gr_hgMapBcType(testBcType,bcType,1,gridDataStruct,axis,face,idest)
     if (testBcType == GRIDBC_GIMME_WORK) then
        call gr_hg_amr_1blk_bcset_work (blockDesc%id, idest, 1, 0)
        return                  ! DONE - RETURN TO CALLER FROM HERE!
     end if
  end if
#endif

  allocate(mask(varCount))
  mask=.true.
  blkLimits   = blockDesc%limits
  blkLimitsGC = blockDesc%limitsGc
  guard=blkLimits(LOW,axis)-blkLimitsGC(LOW,axis)

  isFaceVarNormalDir = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))

  if(regionType(axis)==LEFT_EDGE) then
     sb = blkLimits(LOW,axis)-guard
     se = blkLimits(LOW,axis)+guard-1
     if (isFaceVarNormalDir)se=se+1
  else
     sb = blkLimits(HIGH,axis)-guard+1
     se = blkLimits(HIGH,axis)+guard
     if(isFaceVarNormalDir) sb=sb-1
  end if

  
  if (axis==IAXIS) then
     nextDir(1) = JAXIS
     nextDir(2) = KAXIS
  else if (axis==JAXIS) then
     nextDir(1) = IAXIS
     nextDir(2) = KAXIS
  else
     nextDir(1) = IAXIS
     nextDir(2) = JAXIS
  end if

  do i = 1,2
     j=nextDir(i)
     select case(regionType(j))
     case(LEFT_EDGE)
        endPoints(LOW,j)=blkLimitsGC(LOW,j)
        endPoints(HIGH,j)=blkLimits(LOW,j)-1
     case(CENTER)
        endPoints(:,j)=blkLimits(:,j)
     case(RIGHT_EDGE)
        endPoints(LOW,j)=blkLimits(HIGH,j)+1
        endPoints(HIGH,j)=blkLimitsGC(HIGH,j)
     case(WHOLE_VECTOR)
        endPoints(:,j)=blkLimitsGC(:,j)
     case(NO_VEC)
        endPoints(:,j)=1
     end select
     regionSize(i+1)=endPoints(HIGH,j)-endPoints(LOW,j)+1
  end do

  endPoints(LOW,axis)=sb; endPoints(HIGH,axis)=se
  regionSize(BC_DIR)=se-sb+1
  regionSize(STRUCTSIZE)=varCount
  allocate(regionData(regionSize(BC_DIR),regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),regionSize(STRUCTSIZE)))


  call gr_bcGetRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
       regionData,blockDesc,idest)
  call Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,blockDesc%level,&
       guard,axis,face,regionData,regionSize,mask,applied,&
       nextDir(1),nextDir(2),endPoints,idest)
  if(.not.applied) then
     call Grid_bcApplyToRegion(bcType,gridDataStruct,blockDesc%level,&
          guard,axis,face,regionData,regionSize,mask,applied,&
          nextDir(1),nextDir(2),endPoints,idest)
  end if
  if(.not.applied) then
     print*,'gr_bcApplyToOneFace: Unhandled boundary type',bcType, 'axis,regionType=',axis,regionType(axis)
     if (regionType(axis)==LEFT_EDGE) then
        call Driver_abortFlash("unsupported boundary condition on Lower Face")
     else if (regionType(axis)==RIGHT_EDGE) then
         call Driver_abortFlash("unsupported boundary condition on Upper Face")
     end if
  end if
  call gr_bcPutRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
       regionData,blockDesc,idest)


  deallocate(regionData)
  deallocate(mask)

end subroutine gr_bcApplyToOneFace
