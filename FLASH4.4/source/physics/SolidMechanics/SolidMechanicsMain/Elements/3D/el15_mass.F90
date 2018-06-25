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
  do i=ix,ex
     M(i,i) =  body%mass
  enddo

  ! Orientation variables part:
  select case(trmatrix)

  case( RB_IDENTITY )  

     ! Only iw and ew have values > 0 (see sm_ioRead_rigid.F90), defer.

  case( RB_EULER321 )

     ! NBB and -eye:
     M(ia:ea,ia:ea) = body % RBMAT % NBB
     M(ia,iw)     = -1.
     M(ia+1,iw+1) = -1.
     M(ia+2,iw+2) = -1.

  case( RB_QUATERNN )

     call Driver_abortFlash("el15_mass: error, quaternion mass matrix not coded yet.")

  case default

     call Driver_abortFlash("el15_mass: error, body transformation matrix not known.")

  end select

 
  ! Angular velocity part:
  ! Set Rotational inertia:
  ! Here if TNB=identity -> I_newton = I_body
  TBN = transpose(body%RBMAT%TNB)
  aux_mat = matmul( body%I_body , TBN  )
  body%I_newton = matmul( body%RBMAT%TNB , aux_mat )

  M(iw:ew,iw:ew) = body%I_newton  

  return

end subroutine el15_mass
  
