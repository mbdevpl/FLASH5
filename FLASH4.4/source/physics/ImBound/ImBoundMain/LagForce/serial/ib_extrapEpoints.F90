
subroutine ib_extrapEpoints(nbd,xb,yb,sb,     &
                    ielem,jelem,phile,zL,VarO,&
                    lb,np,lpindex,nx1,ny1,nz1,&
                    del,coord,bsize)


  use ImBound_data , only : ib_ABODY,ib_nmaxa,ib_stencil

  implicit none
#include "Flash.h"
#include "constants.h"
  integer, INTENT(IN) :: nbd,lb,np,nx1,ny1,nz1
  real, INTENT(IN) :: xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa)
  integer, INTENT(IN) :: ielem(ib_stencil,ib_nmaxa),jelem(ib_stencil,ib_nmaxa),lpindex(ib_nmaxa)
  real, INTENT(INOUT) :: phile(ib_stencil,ib_nmaxa)
  real, INTENT(IN) :: zL(ib_nmaxa)
  real, INTENT(INOUT) :: VarO(nx1,ny1,nz1)
  real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)

  ! Local Variables:
  integer :: i,j,ii,ibd,ind1,ind2
  integer :: ipos, iposp1, jpos, jposp1
  real :: dxloc, dyloc, hl, factor
  real :: ds(ib_nmaxa),dsp(ib_nmaxa)
      

  ! Set Force array to zero:
  VarO(:,:,:) = 0.

  ! Define ds, dsp:
  do ibd = 1,nbd

     ind1=ib_ABODY(ibd) % lb + 1; 
     ind2=ib_ABODY(ibd) % lb + ib_ABODY(ibd) % mb; 

     ! Arclength associated with the point
     do ii = ind1,ind2-1
        ds(ii) = sb(ii+1) - sb(ii)
        !ds(ii) = sqrt((xb(ii+1)-xb(ii))**2 + (yb(ii+1)-yb(ii))**2);
     enddo
   
     dsp(ind1) = 0.5*ds(ind1);
     do ii = ind1,ind2-2
        dsp(ii+1) = 0.5*(ds(ii)+ds(ii+1));
     enddo
     dsp(ind2) = 0.5*ds(ind2-1);

  enddo

  ! Local cell size:
  dxloc = del(IAXIS)
  dyloc = del(JAXIS)

  ! Blocks associated Lagrangean points:
  do ii = 1,np

     !factor = 2.*dsp(lpindex(ii))/(dxloc+dyloc); ! Square grids 
     hl = 0.5*(dxloc+dyloc)                       ! Rectangular grids
     factor = dsp(lpindex(ii))*hl/(dxloc*dyloc)

     phile(:,lpindex(ii)) = factor*phile(:,lpindex(ii));
    
     do i = 1,ib_stencil
        
        ipos = ielem(i,lpindex(ii));
        jpos = jelem(i,lpindex(ii));
        
        VarO(ipos,jpos,1) = VarO(ipos,jpos,1) +  &
                            phile(i,lpindex(ii))*zL(lpindex(ii));

     enddo
  enddo

  return

End Subroutine ib_extrapEpoints
