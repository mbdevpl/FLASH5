
#include "ImBound.h"

subroutine ib_getInterpFunc(xp,xyz_stencil,del,derivflag,phile)

  use ImBound_data , only :ib_stencil,ib_interp,ib_npol,ib_alphax,ib_alphay,ib_alphaz

  use ib_interface, only : ib_buildABLan,ib_weightfunc,ib_ludcmp,ib_lubksb

  implicit none
#include "Flash.h"
#include "constants.h"
  real, intent(IN) :: xp(MDIM),xyz_stencil(ib_stencil,MDIM),del(MDIM)
  integer, intent(IN) :: derivflag
  real, intent(OUT):: phile(ib_stencil,NDIM+1)

  ! Local Variables:
  integer, parameter :: nderiv = NDIM+1

  integer :: buildflag

  real,save :: dsx,dsy,dsz

  real    :: A(ib_npol,ib_npol,nderiv), B(ib_npol,ib_stencil,nderiv)
  real    :: p(ib_npol,nderiv), phi(ib_stencil,NDIM+1), gamma(ib_npol,nderiv)
  integer :: indx(ib_npol)
  real    :: d

  logical, save :: firstcall=.true.

  integer :: i,j,idim

  integer :: k,kcount  
  real::xi,yi,zi,wi,wxi,wyi,wzi,wxiyi,wyizi,wxizi,dwidx,dwidy,dwidz
  integer, parameter :: wtype = 1 ! FOR NOW only cubic spline
                                  ! Function W1
  integer INFO

  ! Initialize:
  phile(1:ib_stencil,1:NDIM+1) = 0.

#ifdef FLASH_GRID_PARAMESH
  dsx = ib_alphax*del(IAXIS)
  dsy = ib_alphay*del(JAXIS)
#if NDIM == 3
  dsz = ib_alphaz*del(KAXIS)
#else
  dsz = 0.
#endif
#else
  if (firstcall) then
     dsx = ib_alphax*del(IAXIS)
     dsy = ib_alphay*del(JAXIS)
#if NDIM == 3
     dsz = ib_alphaz*del(KAXIS)
#else
     dsz = 0.
#endif
     firstcall = .false.
  endif
#endif


  ! Build A and B matrices:
#ifndef IB_GET_DERIVS ! NO shape function derivatives are computed in this mode.

  A = 0.
  B = 0.
  buildflag = 0

  k = 0
  kcount = 0
   
  ! Polynomial base p = [1 x y z]
  gamma(1,1) = 1.; gamma(2,1) = xp(IAXIS); gamma(3,1) = xp(JAXIS); gamma(4,1)=xp(KAXIS);

  ! Compute A and B matrices and their derivatives.
  do i = 1 , ib_stencil

      xi=xyz_stencil(i,IAXIS); 
      yi=xyz_stencil(i,JAXIS); 
      zi=xyz_stencil(i,KAXIS);

      !wi,dwidx,dwidy:
      call ib_weightfunc(dsx,dsy,dsz,xp(IAXIS),xp(JAXIS),xp(KAXIS),xi,yi,zi,wtype,wi, &
                         dwidx,dwidy,dwidz)

      ! A and B matrices:
      k = k+1;
      if (wi .NE. 0.) then
         kcount = kcount+1;

         ! Base p =[1 x y z]
         ! A and B:
         wxi  = wi*xi; wyi  = wi*yi; wzi  = wi*zi;
         wxiyi=wxi*yi; wyizi=wyi*zi; wxizi=wxi*zi;

         A(1,1,1) = A(1,1,1) + wi
         A(2,1,1) = A(2,1,1) + wxi
         A(3,1,1) = A(3,1,1) + wyi
         A(4,1,1) = A(4,1,1) + wzi

         A(1,2,1) = A(1,2,1) + wxi
         A(2,2,1) = A(2,2,1) + wxi*xi
         A(3,2,1) = A(3,2,1) + wxiyi
         A(4,2,1) = A(4,2,1) + wxizi

         A(1,3,1) = A(1,3,1) + wyi
         A(2,3,1) = A(2,3,1) + wxiyi
         A(3,3,1) = A(3,3,1) + wyi*yi
         A(4,3,1) = A(4,3,1) + wyizi

         A(1,4,1) = A(1,4,1) + wzi
         A(2,4,1) = A(2,4,1) + wxizi
         A(3,4,1) = A(3,4,1) + wyizi
         A(4,4,1) = A(4,4,1) + wzi*zi

         B(1,k,1) = wi
         B(2,k,1) = wxi
         B(3,k,1) = wyi
         B(4,k,1) = wzi

     endif
  enddo     
  if (kcount .lt. NDIM+1) then !interp=1: in 2D requires kcount 3+, in 3D requires kcount 4+.
     write(*,*) 'kcount',kcount
     A = 0.;
     B = 0.;
     buildflag = 1;
  endif


#else

  call ib_buildABLan(ib_stencil,ib_interp,ib_npol,nderiv,dsx,dsy,dsz,xp(IAXIS),xp(JAXIS),xp(KAXIS),                    &
                     xyz_stencil(1:ib_stencil,IAXIS),xyz_stencil(1:ib_stencil,JAXIS),xyz_stencil(1:ib_stencil,KAXIS),  &
                     A,B,p,buildflag,derivflag); 
  gamma(1:ib_npol,1) = p(1:ib_npol,1)

#endif

  if (buildflag .eq. 1) then
     write(*,*) 'xp, yp, zp=',xp(IAXIS),xp(JAXIS),xp(KAXIS)
     do i=1,ib_stencil
        write(*,*) xyz_stencil(i,IAXIS),xyz_stencil(i,JAXIS),xyz_stencil(i,KAXIS),dsx,dsy,dsz
     enddo
  end if


  ! Obtain gamma coefficients:
  ! Solve for systems coefficients:
  call ib_ludcmp(A(1:ib_npol,1:ib_npol,1),ib_npol,ib_npol,indx,d)

  ! Solve for Gamma:
  !gamma(1:ib_npol,1) = p(1:ib_npol,1)
  call ib_lubksb(A(1:ib_npol,1:ib_npol,1),ib_npol,ib_npol,indx,gamma(1:ib_npol,1))

  ! Obtain Shape functions:
  phi(1:ib_stencil,1) = 0.
  do i = 1 , ib_stencil
     do j = 1, ib_npol
        phi(i,1) = phi(i,1) + gamma(j,1)*B(j,i,1)
     enddo
  enddo
  phile(1:ib_stencil,1) = phi(1:ib_stencil,1);

#ifdef IB_GET_DERIVS
  ! Derivatives of Shape Functions
  if (derivflag .eq. 1) then
     do idim = 1,NDIM

        ! solve for dGamma/dxi:
        gamma(1:ib_npol,idim+1) = p(1:ib_npol,idim+1) - MATMUL(A(1:ib_npol,1:ib_npol,idim+1),gamma(1:ib_npol,1))
        call ib_lubksb(A(1:ib_npol,1:ib_npol,1),ib_npol,ib_npol,indx,gamma(1:ib_npol,idim+1))
        
        ! dphi/dxi:
        phi(1:ib_stencil,idim+1)=MATMUL(gamma(1:ib_npol,idim+1),B(1:ib_npol,1:ib_stencil,1)) + & 
                                 MATMUL(gamma(1:ib_npol,1),B(1:ib_npol,1:ib_stencil,idim+1))

        phile(1:ib_stencil,idim+1) = phi(1:ib_stencil,idim+1)
        
     enddo
  endif
#endif

  return

end subroutine ib_getInterpFunc
