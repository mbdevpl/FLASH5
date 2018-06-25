!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ib_extrapEpoints
!!
!!
!! NAME
!!
!! ib_extrapEpoints
!!
!!
!! SYNOPSIS
!!
!! ib_extrapEpoints(xp,sb,hl,del,ielem,phile,zL,blockID,faceind)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***


subroutine ib_extrapEpoints(xp,sb,hl,del,ielem,phile,zL,blockID,faceind)

  use ImBound_data , only : ib_stencil

  use Grid_interface, ONLY : Grid_getBlkPtr,          &
                             Grid_releaseBlkPtr

  implicit none
#include "Flash.h"
#include "constants.h"
  
  ! Argument list
  integer, INTENT(IN) :: blockID,faceind
  real, INTENT(IN) :: xp(MDIM),sb,hl
  integer, INTENT(IN) :: ielem(ib_stencil,MDIM)
  real, INTENT(IN) :: zL
  real, INTENT(IN) :: del(MDIM)
  real, INTENT(INOUT) :: phile(ib_stencil,NDIM+1)

  ! Local Variables:
  integer :: i
  integer :: ipos, jpos, kpos
  real :: factor
  real, save :: prod
  
  real, pointer, dimension(:,:,:,:) :: faceData

  logical, save :: firstcall=.true.

#ifdef FLASH_GRID_PARAMESH
  prod = PRODUCT(del(1:NDIM))
#else
  if (firstcall) then
     prod = PRODUCT(del(1:NDIM))
     firstcall=.false.
  endif
#endif
  ! Get Pointer to faceData in blockID, direction X,Y,Z:
  call Grid_getBlkPtr(blockID,faceData,faceind)

  ! Rescaling Factor for extrapolation:
  factor = sb*hl/(prod)

  ! Rescale Shape functions:
  ! phile(:,CONSTANT_ONE) = factor*phile(:,CONSTANT_ONE);
    
  ! Do Extrapolation:
  do i = 1,ib_stencil
        
     faceData(FORC_FACE_VAR,ielem(i,IAXIS),ielem(i,JAXIS),ielem(i,KAXIS)) = &
     faceData(FORC_FACE_VAR,ielem(i,IAXIS),ielem(i,JAXIS),ielem(i,KAXIS)) + factor*phile(i,CONSTANT_ONE)*zL;
   
  enddo

  ! Release face data:
  call Grid_releaseBlkPtr(blockID,faceData,faceind)  
  
  return

End Subroutine ib_extrapEpoints
