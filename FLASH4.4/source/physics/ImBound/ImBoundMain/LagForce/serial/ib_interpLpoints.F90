
subroutine ib_interpLpoints(nbd,xb,yb,sb,dsx,dsy, &
         ielem,jelem,phile,zL,VarO,lb,gridflag,np,lpindex, &
         nx1,ny1,nz1,ng,nx,ny,del,coord,bsize)

  use ImBound_data , only : ib_nmaxa,ib_stencil,ib_interp,ib_npol

  use ib_interface , only : ib_buildABLan

  implicit none
#include "Flash.h"
#include "constants.h"

  integer, INTENT(IN) :: nbd,lb,gridflag,np,nx1,ny1,nz1,ng,nx,ny
  real, INTENT(IN) :: xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa)
  real, INTENT(IN) :: dsx(ib_nmaxa),dsy(ib_nmaxa)
  integer, INTENT(IN) :: ielem(ib_stencil,ib_nmaxa),jelem(ib_stencil,ib_nmaxa),lpindex(ib_nmaxa)
  real, INTENT(INOUT) :: phile(ib_stencil,ib_nmaxa),zL(ib_nmaxa)
  real, INTENT(IN) :: VarO(nx1,ny1,nz1)
  real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)

  ! Local Variables:
  integer :: i,j,ii,inod
  integer :: ibd, npoints,buildflag
  real    :: xp,yp,zp,d
  real    :: xi(ib_stencil),yi(ib_stencil)

  real    :: dx,dy,dxaux,dyaux,dVx,dVy
  real    :: A(ib_npol,ib_npol), B(ib_npol,ib_stencil)
  real    :: p(ib_npol), phi(ib_stencil), gamma(ib_npol), indx(ib_npol)

  ! Initialize
  phile(:,:) = 0.0

  dx = del(IAXIS)
  dy = del(JAXIS)

  if(gridflag .eq. IAXIS) then
     dxaux = 0.0
     dyaux = 0.5*dy
  elseif(gridflag .eq. JAXIS) then
     dxaux = 0.5*dx
     dyaux = 0.0
  endif
      

  ! Loop:
  do ibd=1,nbd 

     ! Loop through Eulerian points box:
     do ii = 1,np

        xp = xb(lpindex(ii));
        yp = yb(lpindex(ii));

        xi = coord(IAXIS) - 0.5*bsize(IAXIS) + &
             real(ielem(1:ib_stencil,lpindex(ii)) - ng - 1)*dx + dxaux
        yi = coord(JAXIS) - 0.5*bsize(JAXIS) + &
             real(jelem(1:ib_stencil,lpindex(ii)) - ng - 1)*dy + dyaux
        
        ! Build A and B matrices:
        call ib_buildABLan(ib_stencil,ib_npol,dsx(lpindex(ii)), &
            dsy(lpindex(ii)),xp,yp,xi,yi,ib_interp,A,B,buildflag); 

        ! Obtain gamma coefficients:
        ! Solve for systems coefficients:
        if (ib_interp == 1) then
           p(1) = 1.; p(2) = xp; p(3) = yp;   
        endif

        gamma = p

        if (buildflag .eq. 1) then
           write(*,*) 'xp, yp=',xp,yp
           do i=1,ib_stencil
              write(*,*) xi(i),yi(i),dsx(lpindex(ii)),dsy(lpindex(ii))
           enddo
        end if


        call ludcmp(A,ib_npol,ib_npol,indx,d)
        call lubksb(A,ib_npol,ib_npol,indx,gamma)

        ! Obtain Shape functions:
        phi =0.
        do i = 1 , ib_stencil
           do j = 1, ib_npol
              phi(i) = phi(i) + gamma(j)*B(j,i)
           enddo
        enddo
        phile(1:ib_stencil,lpindex(ii)) = phi(1:ib_stencil);
            
        ! Value of the function in xp and yp:
        zp = 0.;
        do i = 1 , ib_stencil      
           zp = zp + phi(i)* &
                VarO(ielem(i,lpindex(ii)),jelem(i,lpindex(ii)),1);   
        enddo

        zL(lpindex(ii)) = zp;
        
     enddo
  enddo

  return
      
End Subroutine ib_interpLpoints
