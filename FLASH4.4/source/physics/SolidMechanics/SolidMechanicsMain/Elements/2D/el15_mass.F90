!
!
!
! Description : Assemble rigid body Mass matrix
!



subroutine el15_mass(maxdofs,ibd,M)

  use SolidMechanics_Data, only : sm_BodyInfo, sm_structure
  use Driver_interface, only : Driver_abortFlash

  implicit none
#include "SolidMechanics.h"
#include "Flash.h"
#include "constants.h"
  integer, intent(in) :: maxdofs,ibd
  real, intent(out) :: M(maxdofs,maxdofs)

  ! Local variables
  integer :: i,trmatrix
  type(sm_structure), pointer :: body
  real :: TBN(MDIM,MDIM),aux_mat(MDIM,MDIM)
  integer :: ix,ex,ia,ea,iw,ew

  ! Associate Body:
  body => sm_BodyInfo(ibd)

  M   = 0.

  ! Vars to be used:
  ix = body%ix  ! Beginning of translational variables
  ex = body%ex  ! End
  ia = body%ia  ! Beginning of orientation variables
  ea = body%ea  ! End
  iw = body%iw  ! Beginning of angular velocity vars
  ew = body%ew  ! End

  trmatrix = body%trmatrix

  ! M mass matrix related to unrestrained generalized coordinates.
  ! Set values: Translational inertia.
  do i=ix,ex ! In 2D ix=IAXIS, ex=JAXIS
     M(i,i) =  body%mass
  enddo

  ! Orientation variables part:
  select case(trmatrix)

  case( RB_TWODIM )  

     ! Line related ot thetadd: [0 0 1 -1] in 2D ia = 3, ea = 3, iw = 4
     M(ia:ea,ia:ea) =  1.
     M(ia,iw)       = -1.

  case default

     call Driver_abortFlash("el15_mass: error, body transformation matrix not known.")

  end select

 
  ! Angular velocity part:
  ! Set Rotational inertia:
  ! Here -> I_newton = I_body 
  M(iw:ew,iw:ew) = body%I_body(1,1)  

  !write(*,*) M(ix,:)
  !write(*,*) M(ex,:)
  !write(*,*) M(ia,:)
  !write(*,*) M(iw,:)

  return

end subroutine el15_mass
  
