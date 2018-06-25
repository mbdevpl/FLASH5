!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_pk_updatekinematics_rigid
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

subroutine sm_pk_updatekinematics_rigid(ibd,restart,procedure_flag)

  use SolidMechanics_data, only:  sm_BodyInfo, sm_structure
  use sm_pk_data, only: sm_pk_dataset, sm_pk_info
  use sm_pk_interface, only: sm_pk_masterslave_rigid
  use Driver_interface, only : Driver_abortFlash

  implicit none 
#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"
#include "sm_pk.h"
  integer, intent(in) :: ibd,procedure_flag
  logical, intent(in) :: restart

  ! Local Variables:
  type(sm_structure),  pointer :: body
  type(sm_pk_dataset), pointer :: pk
  integer :: nfix, irest, A, i, j, idx
  real, allocatable, dimension(:,:) :: Xi, XiB, v, vd, vdd
  real :: TNB(MDIM,MDIM),NwB_N(MDIM,1),NaB_N(MDIM,1),xp_B(MDIM,1),xp_N(MDIM,1)  
  integer :: maxdofs,nfix_coord,maxrestparams,ircoord
  integer :: imaster
  real :: xm(MDIM,1)
  real, allocatable, dimension(:) :: qn_master,qdn_master,qddn_master


  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Populate orientation angular kinematics vars:
  select case( body % trmatrix )
  case( RB_IDENTITY )
     call sm_NB_kinematicsIdentity(ibd)
     
  case( RB_EULER321 )
     call sm_NB_kinematicsEuler321(ibd)

  case( RB_TWODIM )
     call sm_NB_kinematicsTwoDim(ibd)

  case( RB_QUATERNN )
     call Driver_abortFlash("sm_updatekinematics_rigid: Quaternion technology not coded.")

  case default
     call Driver_abortFlash("sm_updatekinematics_rigid: trmatrix Type not known")
  end select


  ! Initialize or advance:
  if (procedure_flag .eq. SM_INIT) then

     ! Populate initial orientation data for the run:
     select case( body % trmatrix )
     case( RB_IDENTITY , RB_EULER321, RB_TWODIM )

        ! Initial orientation for this run:
        ! Warning : Assume that when simulation is restarted, TNBo
        ! has to be read from the problems initial conditions.
        if (.not. restart) then
           body % RBMAT % TNBo = body % RBMAT % TNB
        else
           body % RBMAT % TNBo = 0.
           if (body % trmatrix .eq. RB_IDENTITY) then
              body % RBMAT % TNBo(1,1) = 1.
              body % RBMAT % TNBo(2,2) = 1.
              body % RBMAT % TNBo(3,3) = 1.
           elseif (body % trmatrix .eq. RB_EULER321) then
              ! Here use problems initial angles.
              body % RBMAT % TNBo(1,1) =  1.
              body % RBMAT % TNBo(2,2) = -1.
              body % RBMAT % TNBo(3,3) = -1.
           elseif (body % trmatrix .eq. RB_TWODIM) then
              ! Here use problems initial angles.
              body % RBMAT % TNBo(1,1) = 1.
              body % RBMAT % TNBo(2,2) = 1.              
           endif
        endif

        ! Fill x,y,z of reference configuration:
        xp_B = 0.
        do i=1,body % nnp
           xp_B(1,1) = body % xB(i)
           xp_B(2,1) = body % yB(i)
#if NDIM == MDIM
           xp_B(3,1) = body % zB(i)
#endif

           ! Position in N frame axes:
           xp_N = matmul(body%RBMAT%TNBo , xp_B) 

           body % x(i) = xp_N(1,1)
           body % y(i) = xp_N(2,1)
#if NDIM == MDIM
           body % z(i) = xp_N(3,1)
#endif

        enddo
     end select

     !do i=1,body % nnp
     !  write(*,*) body % x(i), body % y(i)
     !enddo
     !call Driver_abortFlash("sm_updatekinematics_rigid: After initialize.") 

  elseif(procedure_flag .eq. SM_ADVANCE) then

     ! Now for Master-slave restrained surfaces compute surface node kinematics:
     do irest = 1,body%nrsurf

     ! Get the library of kinematics
     pk => sm_pk_info( body%restraints_surf(irest)%kinematics_idx )

     ! set the number of nodes (not dofs) that are being fixed
     nfix = body%restraints_surf(irest)%nfix

     ! Get the mesh (x,y,z) Xi
     maxdofs   = body%max_dofs_per_node
     allocate( Xi(MDIM,nfix),  XiB(MDIM,nfix), &
               v(maxdofs,nfix), vd(maxdofs,nfix), vdd(maxdofs,nfix) )
     Xi = 0.
     XiB= 0.
     do j = 1,nfix
        A = body%restraints_surf(irest)%node_list(j)
        Xi(1,j) = body%x(A); XiB(1,j) = body%xB(A);
        Xi(2,j) = body%y(A); XiB(2,j) = body%yB(A);
#if NDIM == MDIM
        Xi(3,j) = body%z(A); XiB(3,j) = body%zB(A);
#endif
     end do

     !***************************************
     !* Compute kinematics
     !*
     select case( pk%flag  )

     case( SM_PK_MASTERSLAVE ) ! For rigid body surfaces:
 
        ! Properties of Master Node:
        imaster = int( pk%params(1) )

        ! Extract info from Master:
        ! Reference x,y,z:
        xm = 0.
        xm(1,1) = body % x(imaster)
        xm(2,1) = body % y(imaster)
#if NDIM == MDIM
        xm(3,1) = body % z(imaster)
#endif          

        ! Actual state in qn, qdn, qddn for master Node:
        maxdofs   = body%max_dofs_per_node 
        allocate(qn_master(maxdofs),qdn_master(maxdofs),qddn_master(maxdofs))
        qn_master(1:maxdofs)  = body%  qn(body%ID(1:maxdofs,imaster)) 
        qdn_master(1:maxdofs) = body% qdn(body%ID(1:maxdofs,imaster))
        qddn_master(1:maxdofs)= body%qddn(body%ID(1:maxdofs,imaster))

        ! Orientation variables:
        TNB(1:MDIM,1:MDIM) = body % RBMAT % TNB(1:MDIM,1:MDIM)
        NwB_N(1:MDIM,1)    = body % RBMAT % NwB_N(1:MDIM,1)
        NaB_N(1:MDIM,1)    = body % RBMAT % NaB_N(1:MDIM,1)        

        ! Call Master-Slave constraint imposition on rest Surface: 
        call sm_pk_masterslave_rigid(maxdofs,nfix,xm,qn_master, &
                                     qdn_master,qddn_master,TNB,NwB_N,NaB_N,Xi,XiB,v,vd,vdd)

!!$        i=2
!!$        write(*,*) 'NODE 2 pos=',body%x(i)+body%qn(body%ID(body%ix,i)),   &
!!$                                 body%y(i)+body%qn(body%ID(body%ix+1,i)), &
!!$                                 body%z(i)+body%qn(body%ID(body%ix+2,i))
!!$
!!$        write(*,*) 'NODE 2 vel=',body%qdn(body%ID(body%ix,i)),   &
!!$                                 body%qdn(body%ID(body%ix+1,i)), &
!!$                                 body%qdn(body%ID(body%ix+2,i))
!!$
!!$        write(*,*) 'NODE 2 acc=',body%qddn(body%ID(body%ix,i)),   &
!!$                                 body%qddn(body%ID(body%ix+1,i)), &
!!$                                 body%qddn(body%ID(body%ix+2,i))



        deallocate(qn_master,qdn_master,qddn_master)
 
     end select

     ! Apply kinematics to dofs
     do j = 1,nfix
        A = body%restraints_surf(irest)%node_list(j)
        do i = 1,body%max_dofs_per_node
           idx           = body%ID(i,A)
           body%qn(idx)  = v(i,j)
           body%qdn(idx) = vd(i,j)
           body%qddn(idx)= vdd(i,j)
        end do
     end do

     ! deallocate
     deallocate(Xi, XiB, v, vd, vdd)
     nullify(pk)

     end do

  endif

  return

end subroutine sm_pk_updatekinematics_rigid


!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_NB_kinematicsEuler321
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
subroutine sm_NB_kinematicsEuler321(ibd)

  use SolidMechanics_data, only:  sm_BodyInfo, sm_structure

  implicit none
#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

  integer, intent(in) :: ibd

  ! Local variables:
  type(sm_structure),  pointer :: body
  integer :: maxdofs
  integer :: imaster
  real, allocatable, dimension(:) :: qn_master,qdn_master,qddn_master 
  real :: phi,theta,psi,phid,thetad,psid,phidd,thetadd,psidd

  integer :: ia,iw,ew
  real, parameter :: eps = 1.e-10
  real :: cos_theta
  real, dimension(3,1) :: phiv,phidv,phiddv
  logical :: allpres_ang

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! COM node:
  imaster = Body % Borigin_node

  ! Actual state in qn, qdn, qddn for master Node:
  maxdofs   = body%max_dofs_per_node 
  allocate(qn_master(maxdofs),qdn_master(maxdofs),qddn_master(maxdofs))
  qn_master(1:maxdofs)  = body%  qn(body%ID(1:maxdofs,imaster)) 
  qdn_master(1:maxdofs) = body% qdn(body%ID(1:maxdofs,imaster))
  qddn_master(1:maxdofs)= body%qddn(body%ID(1:maxdofs,imaster))  

  ! phi, theta, psi variables:
  ia = body%ia
  phi   = qn_master(ia);   theta   = qn_master(ia+1);   psi   = qn_master(ia+2);
  phid  = qdn_master(ia);  thetad  = qdn_master(ia+1);  psid  = qdn_master(ia+2);
  phidd = qddn_master(ia); thetadd = qddn_master(ia+1); psidd = qddn_master(ia+2);

  phiv(1,1)   = phi;   phiv(2,1)   = theta;   phiv(3,1)   = psi;
  phidv(1,1)  = phid;  phidv(2,1)  = thetad;  phidv(3,1)  = psid;
  phiddv(1,1) = phidd; phiddv(2,1) = thetadd; phiddv(3,1) = psidd;

  ! This is to avoid singularities on the mass matrix.
  cos_theta = cos(theta)
  if (abs(cos_theta) .le. eps) then
     write(*,*) 'Body=',ibd,'Theta angle close to +-N*pi/2, N odd. cos(theta)=',cos_theta
     cos_theta = eps*sign(1.,cos_theta); 
  end if

  ! Transformation matrix:
  body % RBMAT % TNB         =  0.

  body % RBMAT % TNB(1,1)    =           cos(psi)*cos(theta)
  body % RBMAT % TNB(1,2)    =  sin(phi)*cos(psi)*sin(theta) - cos(phi)*sin(psi)
  body % RBMAT % TNB(1,3)    =  sin(phi)*sin(psi)            + cos(phi)*cos(psi)*sin(theta)  

  body % RBMAT % TNB(2,1)    =          -sin(psi)*cos(theta)
  body % RBMAT % TNB(2,2)    = -sin(phi)*sin(psi)*sin(theta) - cos(phi)*cos(psi)
  body % RBMAT % TNB(2,3)    =  sin(phi)*cos(psi)            - cos(phi)*sin(psi)*sin(theta)  

  body % RBMAT % TNB(3,1)    =                    sin(theta)
  body % RBMAT % TNB(3,2)    = -sin(phi)*         cos(theta)
  body % RBMAT % TNB(3,3)    = -cos(phi)*         cos(theta)


  ! Relation angular velocity - angular dof time derivatives matrices:
  body % RBMAT % NBB         =  0.

  body % RBMAT % NBB(1,1)    =           cos(psi)*cos_theta
  body % RBMAT % NBB(1,2)    = -         sin(psi)
 
  body % RBMAT % NBB(2,1)    = -         sin(psi)*cos_theta
  body % RBMAT % NBB(2,2)    = -         cos(psi)

  body % RBMAT % NBB(3,1)    =                    sin(theta)
  body % RBMAT % NBB(3,3)    = -1.


  body % RBMAT % NBB_dot     =  0.   

  body % RBMAT % NBB_dot(1,1)= -         cos(psi)*sin(theta)*thetad - sin(psi)*cos_theta*psid
  body % RBMAT % NBB_dot(1,2)=                                      - cos(psi)           *psid

  body % RBMAT % NBB_dot(2,1)=           sin(psi)*sin(theta)*thetad - cos(psi)*cos_theta*psid
  body % RBMAT % NBB_dot(2,2)=                                        sin(psi)           *psid

  body % RBMAT % NBB_dot(3,1)=                    cos_theta*thetad  


  ! Angular velocity and acceleration, Newtonian frame N:
  ! There are two options for now:
  ! All prescribed Euler angles:
  iw = body%iw
  ew = body%ew
  allpres_ang = .false.
  allpres_ang = (body%ID(ia,imaster)   .gt. Body%neq) .and. &
                (body%ID(ia+1,imaster) .gt. Body%neq) .and. &
                (body%ID(ia+2,imaster) .gt. Body%neq)


  if (allpres_ang) then

     body % RBMAT % NwB_N = matmul( body%RBMAT%NBB , phidv) 
     body % RBMAT % NaB_N = matmul( body%RBMAT%NBB_dot , phidv) + matmul( body%RBMAT%NBB , phiddv)

     ! Assign to qn,qdn
     body% qdn(body%ID(iw:ew,imaster)) = body % RBMAT % NwB_N(1:3,1)
     body% qddn(body%ID(iw:ew,imaster))= body % RBMAT % NaB_N(1:3,1)


  ! or all free Euler angles
  else
     body % RBMAT % NwB_N(1:3,1) =  qdn_master(iw:ew)
     body % RBMAT % NaB_N(1:3,1) = qddn_master(iw:ew)
  endif

  deallocate(qn_master,qdn_master,qddn_master)

  return

end subroutine sm_NB_kinematicsEuler321

!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_NB_kinematicsIdentity
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
subroutine sm_NB_kinematicsIdentity(ibd)

  use SolidMechanics_data, only:  sm_BodyInfo

  implicit none
#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

  integer, intent(in) :: ibd

  ! Transformation matrix:
  sm_BodyInfo(ibd) % RBMAT % TNB      = 0.
  sm_BodyInfo(ibd) % RBMAT % TNB(1,1) = 1.
  sm_BodyInfo(ibd) % RBMAT % TNB(2,2) = 1.
#if NDIM == MDIM
  sm_BodyInfo(ibd) % RBMAT % TNB(3,3) = 1.
#endif

  ! Angular velocity and acceleration, Newtonian frame N:
  sm_BodyInfo(ibd) % RBMAT % NwB_N = 0.
  sm_BodyInfo(ibd) % RBMAT % NaB_N = 0.

  ! Relation angular velocity - angular dof time derivatives matrices:
  sm_BodyInfo(ibd) % RBMAT % NBB      = 0.
  sm_BodyInfo(ibd) % RBMAT % NBB(1,1) = 1.
  sm_BodyInfo(ibd) % RBMAT % NBB(2,2) = 1.
#if NDIM == MDIM
  sm_BodyInfo(ibd) % RBMAT % NBB(3,3) = 1.
#endif  

  sm_BodyInfo(ibd) % RBMAT % NBB_dot  = 0.   

  return

end subroutine sm_NB_kinematicsIdentity

!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_NB_kinematicsTwoDim
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
subroutine sm_NB_kinematicsTwoDim(ibd)

  use SolidMechanics_data, only:  sm_BodyInfo, sm_structure

  use Driver_interface, only : Driver_abortFlash

  implicit none
#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

  integer, intent(in) :: ibd

  ! Local variables:
  type(sm_structure),  pointer :: body
  integer :: maxdofs
  integer :: imaster
  real, allocatable, dimension(:) :: qn_master,qdn_master,qddn_master 
  real :: theta,thetad,thetadd

  integer :: ia,iw,ew
  real :: cos_theta
  logical :: allpres_ang

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! COM node:
  imaster = Body % Borigin_node

  ! Actual state in qn, qdn, qddn for master Node:
  maxdofs   = body%max_dofs_per_node 
  allocate(qn_master(maxdofs),qdn_master(maxdofs),qddn_master(maxdofs))
  qn_master(1:maxdofs)  = body%  qn(body%ID(1:maxdofs,imaster)) 
  qdn_master(1:maxdofs) = body% qdn(body%ID(1:maxdofs,imaster))
  qddn_master(1:maxdofs)= body%qddn(body%ID(1:maxdofs,imaster))  

  ! Angle theta in 2D:
  ia = body%ia
  theta   = qn_master(ia)  ;   
  thetad  = qdn_master(ia) ; 
  thetadd = qddn_master(ia); 

  ! Transformation matrix:
  body % RBMAT % TNB         =  0.

  body % RBMAT % TNB(1,1)    =           cos(theta)
  body % RBMAT % TNB(1,2)    =          -sin(theta)  

  body % RBMAT % TNB(2,1)    =           sin(theta)
  body % RBMAT % TNB(2,2)    =           cos(theta)


  ! Relation angular velocity - angular dof time derivatives matrices:
  body % RBMAT % NBB         =  0.

  body % RBMAT % NBB(1,1)    =  1.

  body % RBMAT % NBB_dot     =  0.   

  ! Angular velocity and acceleration, Newtonian frame N:
  ! There are two options:
  ! Prescribed orientation angle:
  iw = body%iw
  ew = body%ew
  if ( iw .ne. ew ) &
  call Driver_abortFlash("sm_NB_kinematicsTwoDim: more than one angular velocity scalar in 2D.")
  allpres_ang = .false.
  allpres_ang = (body%ID(ia,imaster)   .gt. Body%neq)

  if (allpres_ang) then

     body % RBMAT % NwB_N(1,1) = thetad;  
     body % RBMAT % NaB_N(1,1) = thetadd; 

     ! Assign to qn,qdn
     body% qdn(body%ID(iw,imaster)) = body % RBMAT % NwB_N(1,1)
     body% qddn(body%ID(iw,imaster))= body % RBMAT % NaB_N(1,1)


  ! or free orientation angle
  else
     body % RBMAT % NwB_N(1,1) =  qdn_master(iw)
     body % RBMAT % NaB_N(1,1) = qddn_master(iw)
  endif

  deallocate(qn_master,qdn_master,qddn_master)

  return

end subroutine sm_NB_kinematicsTwoDim
