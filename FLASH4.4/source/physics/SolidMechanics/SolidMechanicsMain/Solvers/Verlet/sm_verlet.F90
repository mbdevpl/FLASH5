subroutine sm_verlet_partone(ibd,dt)
  
#include "Flash.h"
#include "SolidMechanics.h"
  
  
  use SolidMechanics_rbc_data, ONLY: sm_rbc_Lambda,sm_rbc_mi,n_M
  use SolidMechanics_data, Only : sm_bodyInfo, sm_structure
  
  implicit none
  
  ! Argument list
  integer,INTENT(IN) :: ibd
  real, intent(IN) :: dt
  type(sm_structure),  pointer :: body
  real :: dtt

  body => sm_bodyinfo(ibd)
  ! local arguments
  body%Hs(:) = body%Hs(:)/n_M;
  dtt=0.
  body%qn(:)= body%qms(:,1)  + body%qdms(:,1) * dtt + 0.5 * (body%qddms(:,1) + body%Hs(:)) * dtt *dtt / 1.;
  !body%qn(:)= body%qms(:,1) + body%qdms(:,1) * dt + 
  body%qdi(:)= body%qdms(:,1) + (sm_rbc_Lambda * dtt * (body%qddms(:,1) + body%Hs(:)) ) /1.;

!  Pos2 = Pos1 + (vel1 * dt) + (0.5 * Fn * dt * dt) / 1.; 
!  vel2 = vel1 + (sm_rbc_Lambda * dt * Fn) /1.;
    
  if (maxval(body%qn(:))>100) then
     print*,'There is a problem with the Pos2'
     stop
  end if
end subroutine sm_verlet_partone

subroutine sm_verlet_parttwo(ibd,dt)
  
#include "Flash.h"
#include "SolidMechanics.h"
  
  
  use SolidMechanics_rbc_data, ONLY: sm_rbc_mi, n_M
  use SolidMechanics_data, Only : sm_bodyInfo, sm_structure
  
  implicit none
  
  ! Argument list
  integer,INTENT(IN) :: ibd
  real, intent(IN) :: dt
  type(sm_structure),  pointer :: body
  integer :: i
  real ::ddt
  ddt=0.;
  body => sm_bodyInfo(ibd)
  
  !body%qdms(:,1),body%qdn(:),body%qddms(:,1),body%qddn(:),body%nnp
  body%Hs(:)= body%Hs(:)/n_M
  body%qdn(:)  = body%qdms(:,1) + (0.5 * (body%qddms(:,1)+body%qddn(:)) + body%Hs(:) )*ddt /1.;
  !vel2 = vel1 + (0.5 * dt *(Fn+Fnp1))/1.

  body%qms(:,1)   = body%qn(:);      ! Positions
  body%qdms(:,1)  = body%qdn(:);     ! Velocities
  body%qddms(:,1) = body%qddn(:);    ! Forces


end subroutine sm_verlet_parttwo
