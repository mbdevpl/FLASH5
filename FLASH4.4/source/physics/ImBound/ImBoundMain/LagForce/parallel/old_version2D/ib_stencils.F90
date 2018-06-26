!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/old_version2D/ib_stencils
!!
!!
!! NAME
!!
!! ib_stencils 
!!
!!
!! SYNOPSIS
!!
!!  subroutine ib_stencils(ng,nx,ny,xb,yb,dsx,dsy,   & 
!!                           ielem,jelem,gridflag,  &
!!                           del,coord,bsize)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

subroutine ib_stencils(ng,nx,ny,xb,yb,dsx,dsy,   & 
                       ielem,jelem,gridflag,  &
                       del,coord,bsize)

  use ImBound_data , ONLY : ib_stencil, ib_alphax, ib_alphay
  use gr_sbData, ONLY : gr_sbBodyInfo 

  implicit none
#include "Flash.h"
#include "constants.h"
  
  integer, INTENT(IN) :: ng,nx,ny,gridflag
  real, INTENT(IN) :: xb,yb
  real, INTENT(INOUT) :: dsx,dsy
  integer, INTENT(INOUT) :: ielem(ib_stencil),jelem(ib_stencil)
  real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
   
  ! Local Variables:
  integer :: i,j,k
  real    :: auxx,auxy,distx,disty,xp,yp
  integer :: icpoint,jcpoint
#if NDIM==3
  integer :: kcpoint
#endif
     

  integer :: bndx1,bndx2,bndy1,bndy2,ind1,ind2
  real :: xlower,xupper,ylower,yupper,dx,dy,dxaux,dyaux,xi,yj
  real :: xmin,xmax,ymin,ymax,eps
  real :: aux


  dx = del(IAXIS)
  dy = del(JAXIS)


  eps = 1E-2*MIN(dx,dy)

  bndx1 = ng
  bndx2 = nx     
  bndy1 = ng
  bndy2 = ny     

  dxaux = 0.5*real(gridflag-CONSTANT_ONE)*dx
  dyaux = 0.5*real(CONSTANT_TWO-gridflag)*dy

!!$  if(gridflag .eq. IAXIS) then ! X velocities grid 
!!$     dxaux = 0.0
!!$     dyaux = 0.5*dy
!!$
!!$     ! Find block boundaries
!!$     !xlower = coord(1) - bsize(1)/2.0 - 1.5*dx
!!$     !xupper = coord(1) + bsize(1)/2.0 + 1.5*dx
!!$     !ylower = coord(2) - bsize(2)/2.0 - dy
!!$     !yupper = coord(2) + bsize(2)/2.0 + dy
!!$
!!$  elseif(gridflag .eq. JAXIS) then ! Y velocities grid
!!$     dxaux = 0.5*dx
!!$     dyaux = 0.0
!!$     
!!$     ! Find block boundaries
!!$     !xlower = coord(1) - bsize(1)/2.0 - dx
!!$     !xupper = coord(1) + bsize(1)/2.0 + dx
!!$     !ylower = coord(2) - bsize(2)/2.0 - 1.5*dy
!!$     !yupper = coord(2) + bsize(2)/2.0 + 1.5*dy
!!$     
!!$  endif
  
  
  xp = xb; !
  yp = yb; !

  
  ! Find closest point:
  distx = 10.e10;
  disty = 10.e10;
  icpoint = 0
  jcpoint = 0
  
  ! Closest point in x
  do i = bndx1,bndx2
     xi = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                real(i - ng - 1)*dx + dxaux
     auxx = abs(xp - xi);
     if (auxx+eps .le. distx) then
        icpoint = i;
        distx = auxx;
     endif
  enddo
  
  ! Closest point in y
  do j = bndy1,bndy2
     yj = coord(JAXIS) - 0.5*bsize(JAXIS) + &
          real(j - ng - 1)*dy + dyaux
     auxy = abs(yp - yj);
     if (auxy+eps .le. disty) then
        jcpoint = j;
        disty = auxy;
     endif
  enddo
  
  
  ! dsx and dsy:
  dsx = ib_alphax*dx
  dsy = ib_alphay*dy
  

  ! Build stencil structure:
  ! Center North South East West:
  ielem(1:ib_stencil) = (/ icpoint,icpoint,icpoint,     &
                              icpoint+1,icpoint-1 /)
  jelem(1:ib_stencil) = (/ jcpoint,jcpoint+1,jcpoint-1, &
                              jcpoint,jcpoint /)


  return

end subroutine ib_stencils
