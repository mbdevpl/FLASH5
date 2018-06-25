!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ib_interpLpoints
!!
!!
!! NAME
!!
!! ib_interpLpoints
!!
!!
!! SYNOPSIS
!!
!! ib_interpLpoints(xp,gridfl,del.coord.bsize,   &
!!         ielem,phile,zL,forcflag,blockID,faceind)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

subroutine ib_interpLpoints(xp,gridfl,          &
               del,coord,bsize,ielem,phile,     &
               zL,forcflag,blockID,faceind)

  use ImBound_data , only :ib_stencil,ib_interp,ib_npol

  use ib_interface , only : ib_getInterpFunc

  use Grid_interface, ONLY : Grid_getBlkPtr,          &
                             Grid_releaseBlkPtr

  implicit none
#include "Flash.h"
#include "constants.h"
#include "ImBound.h"

  integer, INTENT(IN) :: gridfl(MDIM),forcflag,blockID,faceind
  real, INTENT(IN)    :: xp(MDIM)
  integer, INTENT(IN) :: ielem(ib_stencil,MDIM)
  real, INTENT(INOUT) :: phile(ib_stencil,NDIM+1),zL
  real, INTENT(IN)    :: del(MDIM),coord(MDIM),bsize(MDIM)

  ! Local Variables:
  integer :: i,idim
  real    :: xyz_stencil(ib_stencil,MDIM)
  real :: delaux(MDIM)
  real :: zp
  integer, parameter :: derivflag = 0 ! Give only Shape functions


  real, pointer, dimension(:,:,:,:) :: faceData


  ! Get Pointer to velocity in blockID, direction X,Y,Z:
  call Grid_getBlkPtr(blockID,faceData,faceind)

  ! Initialize
  delaux(:) = 0.0
  select case (forcflag)
  case(FORCE_FLOW)
     do i=1,NDIM
        delaux(i)   = 0.5*del(i)
        if (gridfl(i) .eq. FACES) delaux(i) = 0.
     enddo
  case(COMPUTE_FORCES)
     do i=1,NDIM
        delaux(i) = 0.
        if (gridfl(i) .eq. CENTER) delaux(i)   = 0.5*del(i)
     enddo
  end select
  

  ! Positions of points on the stencil:
  xyz_stencil(:,:) = 0. 
  do idim = 1,NDIM
     xyz_stencil(:,idim) = coord(idim) - 0.5*bsize(idim) + &
        real(ielem(1:ib_stencil,idim) - NGUARD - 1)*del(idim) + delaux(idim) 
  enddo


  ! Get interpolation functions:
  call ib_getInterpFunc(xp,xyz_stencil,del,derivflag,phile)
         
  ! Value of the function in xp,yp,zp:
  zp = 0.;
  do i = 1 , ib_stencil      
     zp = zp + phile(i,CONSTANT_ONE) * &
               faceData(VELC_FACE_VAR,ielem(i,IAXIS),ielem(i,JAXIS),ielem(i,KAXIS));   
  enddo
  zL = zp;

  ! Release face data (velocities):
  call Grid_releaseBlkPtr(blockID,faceData,faceind)  
  
  return
  
end subroutine ib_interpLpoints
