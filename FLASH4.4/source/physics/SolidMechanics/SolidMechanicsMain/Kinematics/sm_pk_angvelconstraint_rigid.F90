!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_pk_angvelconstraint_rigid
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***

subroutine sm_pk_angvelconstraint_rigid(ibd)

  use SolidMechanics_data, only:  sm_BodyInfo, sm_structure
  use Driver_interface, only : Driver_abortFlash

  implicit none
#include "SolidMechanics.h"
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: ibd
  
  ! Local variables:
  type(sm_structure),  pointer :: body

  integer :: ia,iw,imaster

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Which transformation matrix:
  select case( body % trmatrix )

#if NDIM == MDIM
  case( RB_IDENTITY )

     ! Angular variables all zero:
     ! There is no space for angular vars. Defer.
     
  case( RB_EULER321 )

     call sm_angvelconstraint_Euler321(ibd)

  case( RB_QUATERNN )
     call Driver_abortFlash("sm_pk_angvelconstraint_rigid: Quaternion technology not coded.")

#else /* 2D */

  case( RB_TWODIM )

      ! In 2D thetad and omega are equal by extended mass matrix definition. 
      ! No need to do anything here.
      !ia = body%ia
      !iw = body%iw
      !imaster = Body % Borigin_node 
      !write(*,*) 'Theta=',body%qn(body%ID(ia,imaster)),body%qn(body%ID(iw,imaster))
      !write(*,*) 'Thetad=',body%qdn(body%ID(ia,imaster)),body%qdn(body%ID(iw,imaster))



#endif


  case default
     call Driver_abortFlash("sm_pk_angvelconstraint_rigid: trmatrix Type not known")
  end select  
  
  return

end subroutine sm_pk_angvelconstraint_rigid

!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_angvelconstraint_Euler321
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***
subroutine  sm_angvelconstraint_Euler321(ibd)

  use SolidMechanics_data, only:  sm_BodyInfo, sm_structure

  implicit none
#include "constants.h"
#include "sm_integrator.h"
  integer, intent(in) :: ibd
  
  ! Local variables:
  type(sm_structure),  pointer :: body
  integer :: imaster
  real, parameter :: eps = 1.e-10
  logical :: allfree_ang
  integer :: ia,ea,iw,ew
  real :: theta,psi,cos_theta,iNBB(MDIM,MDIM),NwB_N(1:MDIM,1),phisd(1:MDIM,1)

  ! Get the body and integrator info:
  body  => sm_BodyInfo(ibd)

  ! Only reset orientation time derivatives when the step is converged. Done in sm_PredCorr_advance
  !if( integ%pcconvflag .ne. SM_PCCONVERGED ) return

  ! COM node:
  imaster = Body % Borigin_node

  ! Start and ent of angle and ang vel indexes:
  ia = body % ia
  ea = body % ea
  iw = body % iw
  ew = body % ew

  ! Case all angles are free (moment driven):
  allfree_ang = .false.
  allfree_ang = (body%ID(ia,imaster)   .le. Body%neq) .and. &
                (body%ID(ia+1,imaster) .le. Body%neq) .and. &
                (body%ID(ia+2,imaster) .le. Body%neq)

  if (allfree_ang) then

    theta = body%qn(body%ID(ia+1,imaster)) 
    psi   = body%qn(body%ID(ia+2,imaster))

    ! This is to avoid singularities on inv(NBB).
    cos_theta = cos(theta)
    if (abs(cos_theta) .le. eps) then
      write(*,*) 'Body=',ibd,'Theta angle close to +-N*pi/2, N odd. cos(theta)=',cos_theta
      cos_theta = eps*sign(1.,cos_theta); 
    end if

    ! Compute inverse of NBB:
    iNBB       = 0.
    
    iNBB(1,1)  =  cos(psi)/cos_theta
    iNBB(1,2)  = -sin(psi)/cos_theta
    
    iNBB(2,1)  = -sin(psi)
    iNBB(2,2)  = -cos(psi)
    
    iNBB(3,1)  =  cos(psi)*sin(theta)/cos_theta
    iNBB(3,2)  = -sin(psi)*sin(theta)/cos_theta 
    iNBB(3,3)  = -1.
  
    ! Now get phisd = iNBB * NwB_N
    NwB_N(1:MDIM,1)   = body%qdn(body%ID(iw:ew,imaster));
    phisd   = matmul( iNBB , NwB_N );

    !write(*,*) 'In sm_angvelconstraint_Euler321 B=',body%qdn(body%ID(ia:ea,imaster))

    ! Store in qdn:
    body%qdn(body%ID(ia:ea,imaster))  = phisd(1:MDIM,1);

    !write(*,*) 'In sm_angvelconstraint_Euler321 A=',body%qdn(body%ID(ia:ea,imaster))

  endif

  return

end subroutine sm_angvelconstraint_Euler321
