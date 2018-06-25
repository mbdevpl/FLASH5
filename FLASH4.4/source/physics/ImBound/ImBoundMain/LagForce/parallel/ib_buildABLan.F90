
subroutine ib_buildABLan(stencil,interp,npol,nderiv,dsx,dsy,dsz,xp,yp,zp, &
                         x,y,z,A,B,p,buildflag,derivflag)

  use ib_interface, only : ib_weightfunc

  implicit none
#include "Flash.h"
  integer, INTENT(IN) :: stencil,interp,npol,nderiv,derivflag
  integer, INTENT(OUT) :: buildflag
  real, INTENT(IN) :: dsx,dsy,dsz,xp,yp,zp
  real, INTENT(IN) :: x(stencil),y(stencil),z(stencil)
  real, INTENT(INOUT) :: A(npol,npol,nderiv),B(npol,stencil,nderiv),p(npol,nderiv)

      
  ! Local Variables
  integer :: i,k,kcount
  real :: xi,yi,zi,wi,wxi,wyi,wzi,wxiyi,wyizi,wxizi,dwidx,dwidy,dwidz
  integer, parameter :: wtype = 1 ! FOR NOW only cubic spline
                                  ! Function W1

  buildflag = 0

  A(:,:,:) = 0.
  B(:,:,:) = 0.
  p(:,:) = 0.

  k = 0
  kcount = 0 

  select case (interp)
  case (1)

 
     ! Set Polynomial base:
#if NDIM == 2
     ! Polynomial base p = [1 x y ]
     p(1,1) = 1.; p(2,1) = xp; p(3,1) = yp;
     if (derivflag .eq. 1) then
        ! Polynomial base derivatives:
        ! dp/dx = [0 1 0]
        p(1,2) = 0.; p(2,2) = 1.; p(3,2) = 0.;
        ! dp/dy = [0 0 1]
        p(1,3) = 0.; p(2,3) = 0.; p(3,3) = 1.;
     end if
#elif NDIM == 3
     ! Polynomial base p = [1 x y z]
     p(1,1) = 1.; p(2,1) = xp; p(3,1) = yp; p(4,1)=zp;
     if (derivflag .eq. 1) then
        ! Polynomial base derivatives:
        ! dp/dx = [0 1 0 0]
        p(1,2) = 0.; p(2,2) = 1.; p(3,2) = 0.; p(4,2)=0.;
        ! dp/dy = [0 0 1 0]
        p(1,3) = 0.; p(2,3) = 0.; p(3,3) = 1.; p(4,3)=0.;
        ! dp/dz = [0 0 0 1]
        p(1,4) = 0.; p(2,4) = 0.; p(3,4) = 0.; p(4,4)=1.;
     endif
#endif


     ! Compute A and B matrices and their derivatives.
     do i = 1 , stencil            

        xi=x(i); yi=y(i); zi=z(i)

        !wi,dwidx,dwidy:
        call ib_weightfunc(dsx,dsy,dsz,xp,yp,zp,xi,yi,zi,wtype,wi, &
                           dwidx,dwidy,dwidz)
        
        ! A and B matrices:
        k = k+1;
        if (wi .NE. 0.) then
           kcount = kcount+1;  

#if NDIM == 2 
           ! Base p =[1 x y]
           ! A and B:   
           A(1:3,1,1) = A(1:3,1,1) + wi * &
                      (/ 1.,         x(i),     y(i) /)
           A(1:3,2,1) = A(1:3,2,1) + wi * &
                      (/ x(i),   x(i)**2, x(i)*y(i) /)
           A(1:3,3,1) = A(1:3,3,1) + wi * &
                      (/ y(i), x(i)*y(i),   y(i)**2 /)

           B(1:3,k,1) = (/ wi,x(i)*wi,y(i)*wi/)

#elif NDIM == 3 
           ! Base p =[1 x y z]
           ! A and B:
           wxi  = wi*xi; wyi  = wi*yi; wzi  = wi*zi;
           wxiyi=wxi*yi; wyizi=wyi*zi; wxizi=wxi*zi;
           
           A(1,1,1) = A(1,1,1) + wi
           A(2,1,1) = A(2,1,1) + wxi
           A(3,1,1) = A(3,1,1) + wyi
           A(4,1,1) = A(4,1,1) + wzi
 
!           A(1:4,1,1) = A(1:4,1,1) + wi * &
!                      (/ 1.,         x(i),      y(i),      z(i) /)

           A(1,2,1) = A(1,2,1) + wxi
           A(2,2,1) = A(2,2,1) + wxi*xi
           A(3,2,1) = A(3,2,1) + wxiyi
           A(4,2,1) = A(4,2,1) + wxizi
           
!           A(1:4,2,1) = A(1:4,2,1) + wi * &
!                      (/ x(i),    x(i)**2, x(i)*y(i), x(i)*z(i) /)

           A(1,3,1) = A(1,3,1) + wyi
           A(2,3,1) = A(2,3,1) + wxiyi
           A(3,3,1) = A(3,3,1) + wyi*yi
           A(4,3,1) = A(4,3,1) + wyizi
           
!           A(1:4,3,1) = A(1:4,3,1) + wi * &
!                      (/ y(i), x(i)*y(i),    y(i)**2, y(i)*z(i) /)

           A(1,4,1) = A(1,4,1) + wzi
           A(2,4,1) = A(2,4,1) + wxizi
           A(3,4,1) = A(3,4,1) + wyizi
           A(4,4,1) = A(4,4,1) + wzi*zi

!           A(1:4,4,1) = A(1:4,4,1) + wi * &
!                      (/ z(i), x(i)*z(i),  y(i)*z(i),   z(i)**2 /)
             
           B(1,k,1) = wi
           B(2,k,1) = wxi
           B(3,k,1) = wyi
           B(4,k,1) = wzi

!           B(1:4,k,1) = (/ wi,x(i)*wi,y(i)*wi,z(i)*wi/)

#endif

           if (derivflag .eq. 1) then
              !Derivatives of A and B
#if NDIM == 2

              !dA/dx:
              A(1:3,1,2) = A(1:3,1,2) + dwidx * &
                         (/ 1.,         x(i),     y(i) /)
              A(1:3,2,2) = A(1:3,2,2) + dwidx * &
                         (/ x(i),   x(i)**2, x(i)*y(i) /)
              A(1:3,3,2) = A(1:3,3,2) + dwidx * &
                         (/ y(i), x(i)*y(i),   y(i)**2 /)

              !dA/dy:
              A(1:3,1,3) = A(1:3,1,3) + dwidy * &
                         (/ 1.,         x(i),     y(i) /)
              A(1:3,2,3) = A(1:3,2,3) + dwidy * &
                         (/ x(i),   x(i)**2, x(i)*y(i) /)
              A(1:3,3,3) = A(1:3,3,3) + dwidy * &
                         (/ y(i), x(i)*y(i),   y(i)**2 /)


              !dB/dx, dB/dy:
              B(1:3,k,2) = (/ dwidx,x(i)*dwidx,y(i)*dwidx/)
              B(1:3,k,3) = (/ dwidy,x(i)*dwidy,y(i)*dwidy/)

#elif NDIM == 3

              ! dA/dx:
              A(1:4,1,2) = A(1:4,1,2) + dwidx * &   
                         (/ 1.,         x(i),      y(i),      z(i) /)
              A(1:4,2,2) = A(1:4,2,2) + dwidx * &
                         (/ x(i),    x(i)**2, x(i)*y(i), x(i)*z(i) /)
              A(1:4,3,2) = A(1:4,3,2) + dwidx * &
                         (/ y(i), x(i)*y(i),    y(i)**2, y(i)*z(i) /)
              A(1:4,4,2) = A(1:4,4,2) + dwidx * &
                         (/ z(i), x(i)*z(i),  y(i)*z(i),   z(i)**2 /)
              
              ! dA/dy:
              A(1:4,1,3) = A(1:4,1,3) + dwidy * &   
                         (/ 1.,         x(i),      y(i),      z(i) /)
              A(1:4,2,3) = A(1:4,2,3) + dwidy * &
                         (/ x(i),    x(i)**2, x(i)*y(i), x(i)*z(i) /)
              A(1:4,3,3) = A(1:4,3,3) + dwidy * &
                         (/ y(i), x(i)*y(i),    y(i)**2, y(i)*z(i) /)
              A(1:4,4,3) = A(1:4,4,3) + dwidy * &
                         (/ z(i), x(i)*z(i),  y(i)*z(i),   z(i)**2 /)
             
              ! dA/dz:
              A(1:4,1,4) = A(1:4,1,4) + dwidz * &   
                         (/ 1.,         x(i),      y(i),      z(i) /)
              A(1:4,2,4) = A(1:4,2,4) + dwidz * &
                         (/ x(i),    x(i)**2, x(i)*y(i), x(i)*z(i) /)
              A(1:4,3,4) = A(1:4,3,4) + dwidz * &
                         (/ y(i), x(i)*y(i),    y(i)**2, y(i)*z(i) /)
              A(1:4,4,4) = A(1:4,4,4) + dwidz * &
                         (/ z(i), x(i)*z(i),  y(i)*z(i),   z(i)**2 /)

              ! dB/dx, dB/dy, db/dz:
              B(1:4,k,2) = (/ dwidx,x(i)*dwidx,y(i)*dwidx,z(i)*dwidx/)
              B(1:4,k,3) = (/ dwidy,x(i)*dwidy,y(i)*dwidy,z(i)*dwidy/)
              B(1:4,k,4) = (/ dwidz,x(i)*dwidz,y(i)*dwidz,z(i)*dwidz/)



#endif
           endif
        endif

     enddo

     if (kcount .lt. NDIM+1) then !interp=1: in 2D requires kcount 3+, in 3D requires kcount 4+.
        write(*,*) 'kcount',kcount
        A = 0.;
        B = 0.;
        buildflag = 1;
     endif

  case (2)

     ! Set Polynomial base:
#if NDIM == 2
     ! Polynomial base p = [1 x y x*x x*y y*y]
     p(1:npol,1) = (/ 1., xp, yp, xp**2, xp*yp, yp**2 /)
     if (derivflag .eq. 1) then
        ! Polynomial base derivatives:
        ! dp/dx = [ 0 1 0 2x y 0]
        p(1:npol,2) = (/ 0., 1., 0., 2.*xp, yp, 0. /)
        ! dp/dy = [ 0 0 1 0 x 2y]
        p(1:npol,3) = (/ 0., 0., 1., 0., xp, 2.*yp /) 
     end if
#elif NDIM == 3
     ! Polynomial base p = [1 x y z x*y y*z x*z]
     p(1:npol,1) = (/ 1., xp, yp, zp, xp*yp, yp*zp, xp*zp /)
     if (derivflag .eq. 1) then
        ! Polynomial base derivatives:
        ! dp/dx = [ 0 1 0 0 y 0 z]
        p(1:npol,2) = (/ 0., 1., 0., 0., yp, 0., zp /)
        ! dp/dy = [ 0 0 1 0 x z 0]
        p(1:npol,3) = (/ 0., 0., 1., 0., xp, zp, 0. /)
        ! dp/dz = [ 0 0 0 1 0 y x]
        p(1:npol,4) = (/ 0., 0., 0., 1., 0., yp, xp /)
      endif
#endif

     ! Compute A and B matrices and their derivatives.   
     do i = 1 , stencil

        ! wi,dwidx,dwidy:
        call ib_weightfunc(dsx,dsy,dsz,xp,yp,zp,x(i),y(i),z(i),wtype,wi, &
                           dwidx,dwidy,dwidz);

        
        k = k+1;
        if (wi .NE. 0) then
           kcount = kcount + 1
           ! A and B matrices:
#if NDIM == 2
           ! Polynomial basis p = [1 x y x*x x*y y*y]
           A(1:6,1,1) = A(1:6,1,1) + wi * (/1.,x(i),y(i),x(i)**2,x(i)*y(i),y(i)**2 /)
           A(1:6,2,1) = A(1:6,2,1) + wi * (/x(i),x(i)**2,x(i)*y(i),x(i)**3,x(i)**2*y(i),x(i)*y(i)**2/)
           A(1:6,3,1) = A(1:6,3,1) + wi * (/y(i),x(i)*y(i),y(i)**2,x(i)**2*y(i),x(i)*y(i)**2,y(i)**3/)
           A(1:6,4,1) = A(1:6,4,1) + wi * (/x(i)**2,x(i)**3,x(i)**2*y(i),x(i)**4,x(i)**3*y(i),x(i)**2*y(i)**2/)
           A(1:6,5,1) = A(1:6,5,1) + wi * (/x(i)*y(i),x(i)**2*y(i),x(i)*y(i)**2,x(i)**3*y(i),x(i)**2*y(i)**2,x(i)*y(i)**3/)
           A(1:6,6,1) = A(1:6,6,1) + wi * (/y(i)**2,x(i)*y(i)**2,y(i)**3,x(i)**2*y(i)**2,x(i)*y(i)**3,y(i)**4/)

           B(1:6,k,1) = (/wi,x(i)*wi,y(i)*wi,x(i)**2*wi,x(i)*y(i)*wi,y(i)**2*wi/)

#elif NDIM == 3
           ! Polynomial basis p = [1 x y z x*y y*z x*z] this is not a complete basis.
           A(1:7,1,1) = A(1:7,1,1) + wi*(/ 1., x(i), y(i), z(i), x(i)*y(i),y(i)*z(i),  x(i)*z(i) /)
           A(1:7,2,1) = A(1:7,2,1) + wi*(/x(i),x(i)**2,x(i)*y(i),x(i)*z(i), x(i)**2*y(i), x(i)*y(i)*z(i), x(i)**2*z(i)/)
           A(1:7,3,1) = A(1:7,3,1) + wi*(/y(i),x(i)*y(i),y(i)**2,y(i)*z(i), x(i)*y(i)**2, y(i)**2*z(i), x(i)*y(i)*z(i)/)
           A(1:7,4,1) = A(1:7,4,1) + wi*(/z(i),x(i)*z(i),y(i)*z(i),z(i)**2, x(i)*y(i)*z(i), y(i)*z(i)**2, x(i)*z(i)**2 /)
           A(1:7,5,1) = A(1:7,5,1) + wi*(/ x(i)*y(i), x(i)**2*y(i),x(i)*y(i)**2, x(i)*y(i)*z(i), x(i)**2*y(i)**2, &
                                           x(i)*y(i)**2*z(i), x(i)**2*y(i)*z(i)/)
           A(1:7,6,1) = A(1:7,6,1) + wi*(/ y(i)*z(i), x(i)*y(i)*z(i),y(i)**2*z(i),y(i)*z(i)**2,x(i)*y(i)**2*z(i), &
                                           y(i)**2*z(i)**2,y(i)*z(i)**2*x(i)/)
           A(1:7,7,1) = A(1:7,7,1) + wi*(/ x(i)*z(i), x(i)**2*z(i),x(i)*y(i)*z(i), x(i)*z(i)**2, x(i)**2*y(i)*z(i), &
                                           y(i)*z(i)**2*x(i), x(i)**2*z(i)**2/)
                                
           B(1:7,k,1) = (/wi, x(i)*wi, y(i)*wi, z(i)*wi, x(i)*y(i)*wi,y(i)*z(i)*wi, x(i)*z(i)*wi/)

#endif


           if (derivflag .eq. 1) then
              !Derivatives of A and B
#if NDIM == 2
              !dA/dx:
              A(1:6,1,2) = A(1:6,1,2) + dwidx * (/1.,x(i),y(i),x(i)**2,x(i)*y(i),y(i)**2/)
              A(1:6,2,2) = A(1:6,2,2) + dwidx * (/x(i),x(i)**2,x(i)*y(i),x(i)**3,x(i)**2*y(i),x(i)*y(i)**2/)
              A(1:6,3,2) = A(1:6,3,2) + dwidx * (/y(i),x(i)*y(i),y(i)**2,x(i)**2*y(i),x(i)*y(i)**2,y(i)**3/)
              A(1:6,4,2) = A(1:6,4,2) + dwidx * (/x(i)**2,x(i)**3,x(i)**2*y(i),x(i)**4,x(i)**3*y(i),x(i)**2*y(i)**2/)
              A(1:6,5,2) = A(1:6,5,2) + dwidx * (/x(i)*y(i),x(i)**2*y(i),x(i)*y(i)**2,x(i)**3*y(i),x(i)**2*y(i)**2,x(i)*y(i)**3/)
              A(1:6,6,2) = A(1:6,6,2) + dwidx * (/y(i)**2,x(i)*y(i)**2,y(i)**3,x(i)**2*y(i)**2,x(i)*y(i)**3,y(i)**4/)

              !dA/dy
              A(1:6,1,3) = A(1:6,1,3) + dwidy * (/1.,x(i),y(i),x(i)**2,x(i)*y(i),y(i)**2 /)
              A(1:6,2,3) = A(1:6,2,3) + dwidy * (/x(i),x(i)**2,x(i)*y(i),x(i)**3,x(i)**2*y(i),x(i)*y(i)**2/)
              A(1:6,3,3) = A(1:6,3,3) + dwidy * (/y(i),x(i)*y(i),y(i)**2,x(i)**2*y(i),x(i)*y(i)**2,y(i)**3/)
              A(1:6,4,3) = A(1:6,4,3) + dwidy * (/x(i)**2,x(i)**3,x(i)**2*y(i),x(i)**4,x(i)**3*y(i),x(i)**2*y(i)**2/)
              A(1:6,5,3) = A(1:6,5,3) + dwidy * (/x(i)*y(i),x(i)**2*y(i),x(i)*y(i)**2,x(i)**3*y(i),x(i)**2*y(i)**2,x(i)*y(i)**3/)
              A(1:6,6,3) = A(1:6,6,3) + dwidy * (/y(i)**2,x(i)*y(i)**2,y(i)**3,x(i)**2*y(i)**2,x(i)*y(i)**3,y(i)**4/)

              !dB/dx and dB/dy
              B(1:6,k,2) = (/dwidx,x(i)*dwidx,y(i)*dwidx,x(i)**2*dwidx,x(i)*y(i)*dwidx,y(i)**2*dwidx/)
              B(1:6,k,3) = (/dwidy,x(i)*dwidy,y(i)*dwidy,x(i)**2*dwidy,x(i)*y(i)*dwidy,y(i)**2*dwidy/)

#elif NDIM == 3
              ! dA/dx:
              A(1:7,1,2) = A(1:7,1,2)+dwidx*(/ 1., x(i), y(i), z(i),x(i)*y(i), &
                                             y(i)*z(i),  x(i)*z(i) /)

              A(1:7,2,2) = A(1:7,2,2)+dwidx*(/x(i),x(i)**2,x(i)*y(i),x(i)*z(i), &
                                             x(i)**2*y(i), x(i)*y(i)*z(i), x(i)**2*z(i)/)

              A(1:7,3,2) = A(1:7,3,2)+dwidx*(/y(i),x(i)*y(i),y(i)**2,y(i)*z(i), &
                                             x(i)*y(i)**2, y(i)**2*z(i), x(i)*y(i)*z(i)/)

              A(1:7,4,2) = A(1:7,4,2)+dwidx*(/z(i),x(i)*z(i),y(i)*z(i),z(i)**2, &
                                             x(i)*y(i)*z(i), y(i)*z(i)**2, x(i)*z(i)**2 /)

              A(1:7,5,2) = A(1:7,5,2)+dwidx*(/ x(i)*y(i), x(i)**2*y(i), &
                                             x(i)*y(i)**2, x(i)*y(i)*z(i), x(i)**2*y(i)**2, &
                                             x(i)*y(i)**2*z(i), x(i)**2*y(i)*z(i)/)

              A(1:7,6,2) = A(1:7,6,2)+dwidx*(/ y(i)*z(i), x(i)*y(i)*z(i),y(i)**2*z(i), & 
                                             y(i)*z(i)**2,x(i)*y(i)**2*z(i),y(i)**2*z(i)**2, &
                                             y(i)*z(i)**2*x(i)/)

              A(1:7,7,2) = A(1:7,7,2)+dwidx*(/ x(i)*z(i), x(i)**2*z(i), x(i)*y(i)*z(i), &
                                             x(i)*z(i)**2, x(i)**2*y(i)*z(i), y(i)*z(i)**2*x(i), &
                                             x(i)**2*z(i)**2/)


              ! dA/dy: 
              A(1:7,1,3) = A(1:7,1,3)+dwidy*(/ 1., x(i), y(i), z(i),x(i)*y(i), &
                                             y(i)*z(i),  x(i)*z(i) /)

              A(1:7,2,3) = A(1:7,2,3)+dwidy*(/x(i),x(i)**2,x(i)*y(i),x(i)*z(i), &
                                             x(i)**2*y(i), x(i)*y(i)*z(i), x(i)**2*z(i)/)

              A(1:7,3,3) = A(1:7,3,3)+dwidy*(/y(i),x(i)*y(i),y(i)**2,y(i)*z(i), &
                                             x(i)*y(i)**2, y(i)**2*z(i), x(i)*y(i)*z(i)/)

              A(1:7,4,3) = A(1:7,4,3)+dwidy*(/z(i),x(i)*z(i),y(i)*z(i),z(i)**2, &
                                             x(i)*y(i)*z(i), y(i)*z(i)**2, x(i)*z(i)**2 /)

              A(1:7,5,3) = A(1:7,5,3)+dwidy*(/ x(i)*y(i), x(i)**2*y(i), &
                                             x(i)*y(i)**2, x(i)*y(i)*z(i), x(i)**2*y(i)**2, &
                                             x(i)*y(i)**2*z(i), x(i)**2*y(i)*z(i)/)

              A(1:7,6,3) = A(1:7,6,3)+dwidy*(/ y(i)*z(i), x(i)*y(i)*z(i),y(i)**2*z(i), & 
                                             y(i)*z(i)**2,x(i)*y(i)**2*z(i),y(i)**2*z(i)**2, &
                                             y(i)*z(i)**2*x(i)/)

              A(1:7,7,3) = A(1:7,7,3)+dwidy*(/ x(i)*z(i), x(i)**2*z(i), x(i)*y(i)*z(i), &
                                             x(i)*z(i)**2, x(i)**2*y(i)*z(i), y(i)*z(i)**2*x(i), &
                                             x(i)**2*z(i)**2/)


            
              ! dA/dz:
              A(1:7,1,4) = A(1:7,1,4)+dwidz*(/ 1., x(i), y(i), z(i),x(i)*y(i), &
                                             y(i)*z(i),  x(i)*z(i) /)

              A(1:7,2,4) = A(1:7,2,4)+dwidz*(/x(i),x(i)**2,x(i)*y(i),x(i)*z(i), &
                                             x(i)**2*y(i), x(i)*y(i)*z(i), x(i)**2*z(i)/)

              A(1:7,3,4) = A(1:7,3,4)+dwidz*(/y(i),x(i)*y(i),y(i)**2,y(i)*z(i), &
                                             x(i)*y(i)**2, y(i)**2*z(i), x(i)*y(i)*z(i)/)

              A(1:7,4,4) = A(1:7,4,4)+dwidz*(/z(i),x(i)*z(i),y(i)*z(i),z(i)**2, &
                                             x(i)*y(i)*z(i), y(i)*z(i)**2, x(i)*z(i)**2 /)

              A(1:7,5,4) = A(1:7,5,4)+dwidz*(/ x(i)*y(i), x(i)**2*y(i), &
                                             x(i)*y(i)**2, x(i)*y(i)*z(i), x(i)**2*y(i)**2, &
                                             x(i)*y(i)**2*z(i), x(i)**2*y(i)*z(i)/)

              A(1:7,6,4) = A(1:7,6,4)+dwidz*(/ y(i)*z(i), x(i)*y(i)*z(i),y(i)**2*z(i), & 
                                             y(i)*z(i)**2,x(i)*y(i)**2*z(i),y(i)**2*z(i)**2, &
                                             y(i)*z(i)**2*x(i)/)

              A(1:7,7,4) = A(1:7,7,4)+dwidz*(/ x(i)*z(i), x(i)**2*z(i), x(i)*y(i)*z(i), &
                                             x(i)*z(i)**2, x(i)**2*y(i)*z(i), y(i)*z(i)**2*x(i), &
                                             x(i)**2*z(i)**2/)



              ! dB/dx, dB/dy, db/dz:
              B(1:7,k,2) = (/ dwidx,x(i)*dwidx,y(i)*dwidx,z(i)*dwidx,x(i)*y(i)*dwidx, y(i)*z(i)*dwidx, x(i)*z(i)*dwidx/)
              B(1:7,k,3) = (/ dwidy,x(i)*dwidy,y(i)*dwidy,z(i)*dwidy,x(i)*y(i)*dwidy, y(i)*z(i)*dwidy, x(i)*z(i)*dwidy/)
              B(1:7,k,4) = (/ dwidz,x(i)*dwidz,y(i)*dwidz,z(i)*dwidz,x(i)*y(i)*dwidz, y(i)*z(i)*dwidz, x(i)*z(i)*dwidz/)

#endif

           endif

        endif
     enddo

     if (kcount .lt. NDIM+4) then !interp=2: in 2D requires kcount 6+, in 3D requires kcount 7+. 
        A = 0.;
        B = 0.;
        buildflag = 1;
     endif

  end select

  return

End subroutine ib_buildABLan
