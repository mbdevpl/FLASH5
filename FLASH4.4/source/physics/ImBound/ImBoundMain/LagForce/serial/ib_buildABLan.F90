
subroutine ib_buildABLan(stencil,npol,dsx,dsy,xp,yp, &
                         x,y,interp,A,B,buildflag)

  use ib_interface, only : ib_weightfunc

  implicit none
  integer, INTENT(IN) :: stencil,npol,interp
  integer, INTENT(OUT) :: buildflag
  real, INTENT(IN) :: dsx,dsy,xp,yp
  real, INTENT(IN) :: x(stencil),y(stencil)
  real, INTENT(INOUT) :: A(npol,npol),B(npol,stencil)

      
  ! Local Variables
  integer :: i,k,kcount
  real :: wi,dwidx,dwidy
  integer, parameter :: wtype = 1 ! FOR NOW only cubic spline
                                  ! Function W1

  buildflag = 0

  select case (interp)
  case (1)

     A = 0.;
     B = 0.;

     k = 0;
     kcount = 0; 
        
     do i = 1 , stencil            

        !wi,dwidx,dwidy:
        call ib_weightfunc(dsx,dsy,xp,yp,x(i),y(i),wtype,wi, &
                           dwidx,dwidy)
        
        ! A and B matrices:
        k = k+1;
        if (wi .NE. 0.) then
           kcount = kcount +1;                
           A(1:3,1) = A(1:3,1) + wi * &
                (/ 1.,         x(i),      y(i) /)
           A(1:3,2) = A(1:3,2) + wi * &
                (/ x(i),   x(i)**2, x(i)*y(i) /)
           A(1:3,3) = A(1:3,3) + wi * &
                (/ y(i), x(i)*y(i),   y(i)**2 /)

           B(:,k) = (/ wi,x(i)*wi,y(i)*wi/)
        endif

     enddo

     if (kcount < 3) then
        write(*,*) 'kcount',kcount
        A = 0.;
        B = 0.;
        buildflag = 1;
     endif

  case (2)

     A = 0.;
     B = 0.;

     k = 0;
     kcount = 0;
        
     do i = 1 , stencil

        ! wi,dwidx,dwidy:
        call ib_weightfunc(dsx,dsy,xp,yp,x(i),y(i),wtype,wi, &
                           dwidx,dwidy);

        ! A and B matrices:
        k = k+1;
        if (wi .NE. 0) then
           kcount = kcount + 1
           A = A + wi * reshape((/1.,x(i),y(i),x(i)**2,          &
                   x(i)*y(i),y(i)**2,x(i),x(i)**2,               &
                   x(i)*y(i),x(i)**3,x(i)**2*y(i),x(i)*y(i)**2,  &
                   y(i),x(i)*y(i),y(i)**2,x(i)**2*y(i),          &
                   x(i)*y(i)**2,y(i)**3,x(i)**2,x(i)**3,         &
                   x(i)**2*y(i),x(i)**4,x(i)**3*y(i),            &
                   x(i)**2*y(i)**2,x(i)*y(i),x(i)**2*y(i),       &
                   x(i)*y(i)**2,x(i)**3*y(i),x(i)**2*y(i)**2,    &
                   x(i)*y(i)**3,y(i)**2,x(i)*y(i)**2,y(i)**3,    &
                   x(i)**2*y(i)**2,x(i)*y(i)**3,y(i)**4/),       &
                   (/6,6/))

           B(:,k) = (/wi,x(i)*wi,y(i)*wi,x(i)**2*wi, &
                         x(i)*y(i)*wi,y(i)**2*wi/)

        endif
     enddo

     if (kcount < 6) then
        A = 0.;
        B = 0.;
        buildflag = 1;
     endif

  end select

  return

End subroutine ib_buildABLan
