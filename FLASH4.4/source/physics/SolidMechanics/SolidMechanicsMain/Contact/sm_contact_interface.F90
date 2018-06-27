module sm_contact_interface
  interface	
     function cross(a,b)
       real ,dimension(1,3)::a,b,cross
     end function cross
  end interface
  interface 
     function length(a)
       real ,dimension(1,3)::a
       real :: length
     end function length
  end interface
  interface 
     subroutine  sm_contact(ibd,restart_local)
       
#include "SolidMechanics.h"
#include "Flash.h"
      
       implicit none
       integer, intent(IN) :: ibd
       logical,intent(IN):: restart_local
     end subroutine sm_contact
  end interface

  interface 
     subroutine collisionVolume(xyz2,ind,intVol)
       
#include "constants.h"
       
       implicit none
       integer, intent(IN) :: ind
       real, intent(IN),dimension(ind,3):: xyz2
       real, intent(OUT) :: intVol
     end subroutine collisionVolume
  end interface

  interface
     subroutine detectCollision(point,vel,s1,s2,s3,s_vel1,s_vel2,s_vel3,intersect,dt_local)
       implicit none
       real, dimension(1,3), intent(IN)::point,vel,s1,s2,s3,s_vel1,s_vel2,s_vel3
       real, intent(IN) :: dt_local
       real, intent(INOUT) :: intersect
       
     end subroutine detectCollision
  end interface


  interface
     subroutine sm_contact_init()
       implicit none
     end subroutine sm_contact_init
  end interface


end module sm_contact_interface


function length(a) 
  real ,dimension(3)::a
  real ::length

  length=SQRT(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
  
  return
end function length




function cross(a,b)
  real ,dimension(1,3)::a,b,cross
  
  cross(1,1)=a(1,2)*b(1,3)-a(1,3)*b(1,2)
  cross(1,2)=a(1,3)*b(1,1)-a(1,1)*b(1,3)
  cross(1,3)=a(1,1)*b(1,2)-a(1,2)*b(1,1)
  
  return
end function cross

function isnan(a)
  real ,dimension(1,1)::a
  logical::isnan
  
  if(a(1,1) .ne. a(1,1)) then
     isnan = .true.
  else
     isnan = .false.
  endif
  return
end function isnan

function isnan2(a)
  real::a
  logical::isnan2
  
  if(a .ne. a) then
     isnan2= .true.
  else
     isnan2= .false.
  endif
  return
end function isnan2
