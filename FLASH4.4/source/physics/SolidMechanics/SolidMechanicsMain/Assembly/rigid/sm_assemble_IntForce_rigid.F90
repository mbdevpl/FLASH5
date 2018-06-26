!    
! File:  
! 
! Inertia driven RHS for Rigid EOMS.

subroutine sm_assemble_IntForce_rigid(ibd, flag)

  use SolidMechanics_Data, only : sm_BodyInfo, sm_structure
  use sm_Misc_interface, only : sm_crossProd
  use Driver_interface, only : Driver_abortFlash

  implicit none
#include "SolidMechanics.h"
#include "Flash.h"
#include "constants.h"    
  !IO Variables
  integer, intent(in)    :: ibd, flag ! body number
  
  ! Local variables:
  integer :: trmatrix
  type(sm_structure), pointer :: body
  real :: TBN(MDIM,MDIM),aux_mat(MDIM,MDIM),aux_v(MDIM,1),phisd(MDIM,1)
  integer :: ix,ex,ia,ea,iw,ew,imaster,idx,i,i_dim
  
  real, allocatable, dimension(:) :: Qint
  real, allocatable, dimension(:) :: qn_master


  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Vars to be used:
  ix = body%ix  ! Beginning of translational variables
  ex = body%ex  ! End
  ia = body%ia  ! Beginning of orientation variables
  ea = body%ea  ! End
  iw = body%iw  ! Beginning of angular velocity vars
  ew = body%ew  ! End

  trmatrix = body%trmatrix

  ! Master node:
  imaster = body%borigin_node

  ! Initialize:
  body%Qs = 0.
  allocate(Qint(body%max_dofs_per_node)); Qint=0.
  allocate(qn_master(body%max_dofs_per_node)); qn_master=0.

#if NDIM==MDIM

  ! Set Rotational inertia:
  ! Here if TNB=identity -> I_newton = I_body
  TBN = transpose(body%RBMAT%TNB)
  aux_mat = matmul( body%I_body , TBN  )
  body%I_newton = matmul( body%RBMAT%TNB , aux_mat )

  ! Internal forces - RB inertia driven:
  ! Add rotation part:
  select case(trmatrix)

  case( RB_IDENTITY )  

     !  defer.

  case( RB_EULER321 )

     ! -NBB_dot * phisd:
     phisd(1,1) = body%qdn(body%ID(ia  ,imaster))
     phisd(2,1) = body%qdn(body%ID(ia+1,imaster))
     phisd(3,1) = body%qdn(body%ID(ia+2,imaster))

     aux_v = -matmul( body%RBMAT%NBB_dot , phisd)

     ! Add to Qint.
     i_dim = 0
     do i = ia,ea
        i_dim = i_dim + 1
        Qint(i) = Qint(i) + aux_v(i_dim,1)
     end do

  case( RB_QUATERNN )

     call Driver_abortFlash("sm_assemble_Intforce_rigid: error, quaternion vector not coded yet.")

  case default

     call Driver_abortFlash("sm_assemble_Intforce_rigid: error, body transformation matrix not known.")

  end select  

  ! Add Angular velocity part - NwB x IN*NwB :
  aux_v = matmul( body%I_newton , body%RBMAT%NwB_N )
  Qint(iw:ew) = - sm_crossProd( body%RBMAT%NwB_N(1:NDIM,1) , aux_v(1:NDIM,1) )

#else /* 2D */

  ! Case 2D add stiffness vector:
  qn_master(ix:ew) = body%qn(body%ID(ix:ew,imaster))

  Qint(ix)    = - body%Stiff(IAXIS) * qn_master(ix)
  Qint(ex)    = - body%Stiff(JAXIS) * qn_master(ex)
  Qint(ia)    = 0.  ! ea = ia = 3 in 2D
  Qint(iw)    = - body%Stiff(ia) * qn_master(ia) ! wdot line -Ktheta*theta, iw = ew in 2D

#endif

  ! Add to whatever is in Qs:
  ! Forces:
  do i = body%ix,body%ex
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Qs(idx)  = body%Qs(idx) + Qint(i)
     endif
  end do

  ! Orientation:
  do i = body%ia,body%ea
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Qs(idx)  = body%Qs(idx) + Qint(i)
     endif
  end do 

  ! Moments: 
  do i = body%iw,body%ew
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Qs(idx)  = body%Qs(idx) + Qint(i)
     endif
  end do

  deallocate(Qint,qn_master)

  return
    
end subroutine sm_assemble_IntForce_rigid
