subroutine el04_IntForce_rbc (p1,p2,p3,p4,  &
     n1,n2,c1,c2,a1,a2,th,        &
     tho,cos_tho,sin_tho, bforces)


#include "Flash.h"
#include "SolidMechanics.h"

  use SolidMechanics_rbc_data, only : kbf_M
  
  implicit none
  interface
     function cross_product(a,b)
       real ::a(3),b(3),cross_product(3)
     end function cross_product
  end interface
  
  
  real,dimension(NDIM),intent(IN)::p1,p2,p3,p4
  real,dimension(NDIM),intent(IN)::n1,n2,c1,c2
  real, intent(IN) :: a1, a2
  real, intent(IN) :: tho,cos_tho,sin_tho
  real, dimension(MAXNODERBC,NDIM), intent(OUT) :: bforces
  real, intent(OUT) :: th
  
  real,dimension(NDIM)::a32,a13,a34,a21,a42,a23
  real::b11,b12,b22,beta_b,nm1,nm2     
  real::coeff_sign
  real,dimension(NDIM)::n1a32,n1a13,n1a34,n1a42
  real,dimension(NDIM)::n1a23,n2a32,n2a13,n2a34
  real,dimension(NDIM)::n2a42,n2a23,n1a21,n2a21  
  real:: sinsign, sin_th, cos_th
 
  
  a32=p3-p2
  a13=p1-p3
  a34=p3-p4
  a21=p2-p1
  a42=p4-p2
  a23=p2-p3
  
  n1a32=cross_product(n1,a32)
  n1a13=cross_product(n1,a13)
  n1a34=cross_product(n1,a34)
  n1a21=cross_product(n1,a21)
  n1a42=cross_product(n1,a42)
  n1a23=cross_product(n1,a23)
  
  n2a32=cross_product(n2,a32)
  n2a13=cross_product(n2,a13)
  n2a34=cross_product(n2,a34)
  n2a21=cross_product(n2,a21)
  n2a42=cross_product(n2,a42)
  n2a23=cross_product(n2,a23)
  
  nm1=a1*2.
  nm2=a2*2.
 
  cos_th=dot_product((n1/nm1),(n2/nm2))
  
  
  if (cos_th.gt.0.999999500000042) then
     
     cos_th=0.999999500000042
     
  elseif(cos_th.LE.-1.0) then
   
     cos_th=-1.0;
     write(*,*) cos_th
   
     write(*,*) 'p1,p2',p1,p2
     write(*,*) 'p3,p4', p3,p4
     write(*,*) 'n1,n2',n1,n2
     write(*,*) 'c1,c2',c1,c2
     write(*,*) 'a1,a2,th,',a1,a2,th

  endif
  
  sin_th=sqrt(1.-(cos_th*cos_th));
  th=acos(cos_th)
  !write(*,*) 'th',th

  sinsign=dot_product((n1-n2),(c1-c2))
  IF (sinsign.LE.0) sin_th=-1.*sin_th;
  !write(*,*) kbf_M, cos_tho, cos_th, sin_tho, sin_th
  
  beta_b=kbf_M*(sin_th*cos_tho-cos_th*sin_tho)/sin_th
  
  b11=-1.*beta_b*cos_th/(nm1*nm1)
  
  b12=beta_b/(nm1*nm2)
  
  b22=-1.*beta_b*cos_th/(nm2*nm2)
  
  bforces(1,1:3)=( b11*n1a32 +b12* n2a32)
  bforces(2,1:3)=((b11*n1a13)+b12*(n1a34+n2a13)+(b22*n2a34))
  bforces(3,1:3)=((b11*n1a21)+b12*(n1a42+n2a21)+(b22*n2a42))
  bforces(4,1:3)=( b12*n1a23 +b22* n2a23);
  
  if (maxval(abs(bforces(:,:)))>1000) then
     write(*,*) 'beta_b=',beta_b,b11,b12,b22
     write(*,*) ' nm1,nm2',nm1,nm2
     print*, 'Problem in the bending forces'
  !   stop
  end if
  
end subroutine el04_IntForce_rbc

