subroutine el03_IntForce_rbc(p1,p2,p3,  & 
     area,Aoj,At, Aot, &
     Vt, Vot, normal,center, tforces)
  
#include "Flash.h"
#include "SolidMechanics.h"
  
  use SolidMechanics_rbc_data, only: sm_rbc_ka,sm_rbc_kv,sm_rbc_kd 
  
  implicit none
  
  Interface
     function cross_product(a,b)   
       real ::a(3),b(3),cross_product(3)        
     end function cross_product
  end interface
  
  real,dimension(NDIM), intent(IN) :: p1,p2,p3
  real,dimension(NDIM), intent(IN) :: normal,center
  real,dimension(MAXNODERBC,NDIM),intent(OUT)::tforces
  real,intent(IN) ::Aoj,Vt,area,At
  real,intent(IN) ::Aot,Vot
  
  ! Local variables
  real ::beta_a,beta_v,beta_d
  real,dimension(NDIM)::a32,a13,a21
  real::alpha
  real,dimension(NDIM)::func1,func2,func3
  real,dimension(9,NDIM)::ttforces

  a32=p3-p2
  a13=p1-p3
  a21=p2-p1
  ttforces=0.0;

  ! functional form
  func1=cross_product(normal,a32)
  func2=cross_product(normal,a13)
  func3=cross_product(normal,a21)
!!$  write(*,*)func1
!!$  write(*,*)func2
!!$  write(*,*)func3
  ! Global area effect
  alpha=-1.*real(sm_rbc_ka)*(At-Aot)/(Aot*4.*area)
  
  ttforces(1,1:3)=alpha*func1
  ttforces(2,1:3)=alpha*func2
  ttforces(3,1:3)=alpha*func3
  
  ! local area constraint forces
  beta_d=-1.*real(sm_rbc_kd)*(area-Aoj)/(Aoj*4.*area)
  

  ttforces(4,1:3)=beta_d*func1
  ttforces(5,1:3)=beta_d*func2
  ttforces(6,1:3)=beta_d*func3

  !     volume constraint forces
  beta_v=-1.*real(sm_rbc_kv)*(Vt-Vot)/(Vot*6.)
!!$  write(*,*)alpha,beta_d,beta_v
!!$  write(*,*)normal
!!$  write(*,*)center
!!$  write(*,*)cross_product(center,a32)
!!$  write(*,*)cross_product(center,a13)
!!$  write(*,*)cross_product(center,a21)
    

  ttforces(7,1:3)=beta_v*(normal/3.+cross_product(center,a32))
  ttforces(8,1:3)=beta_v*(normal/3.+cross_product(center,a13))
  ttforces(9,1:3)=beta_v*(normal/3.+cross_product(center,a21))

  
  tforces(1,1:3)= ttforces(1,1:3) + ttforces(4,1:3) + ttforces(7,1:3)
  tforces(2,1:3)= ttforces(2,1:3) + ttforces(5,1:3) + ttforces(8,1:3)
  tforces(3,1:3)= ttforces(3,1:3) + ttforces(6,1:3) + ttforces(9,1:3)
!!$  write(*,*)tforces(1,1:3)
!!$  write(*,*)tforces(2,1:3)
!!$  write(*,*)tforces(3,1:3)
  !stop

  if (maxval(abs(tforces(1:3,1:3)))>2000) then
     write(*,*)p1
     write(*,*)p2
     write(*,*)p3
     write(*,*)'norm',normal
     write(*,*)'cent',center
     write(*,*) 'Aoj',Aoj,'Area',area
     write(*,*)ttforces(1,1:3)
     write(*,*)ttforces(2,1:3)
     write(*,*)ttforces(3,1:3)
     write(*,*)ttforces(4,1:3)
     write(*,*)ttforces(5,1:3)
     write(*,*)ttforces(6,1:3)
     write(*,*)ttforces(7,1:3)
     write(*,*)ttforces(8,1:3)
     write(*,*)ttforces(9,1:3)
     write(*,*) 'qes'
     write(*,*)tforces(1,1:3)
     write(*,*)tforces(2,1:3)
     write(*,*)tforces(3,1:3)
     print*, 'Problem in the area forces'
  !   stop
  end if
end subroutine el03_IntForce_rbc

