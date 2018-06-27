!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_pk_apply_entireBody
!!
!!
!! NAME
!!
!! 
!!
!!
!! SYNOPSIS
!!    apply the prescribed kinematics of RestSurface 1 to the entire body
!!  
!!
!!
!! DESCRIPTION
!!  
!!
!!
!!***

#include "sm_pk.h"
#include "Flash.h"
#include "constants.h"

subroutine sm_pk_apply_entireBody(ibd,time)

  use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
  use sm_pk_data, only: sm_pk_dataset, sm_pk_info
  use sm_pk_interface, only: sm_pk_Flapping_BermanWang, sm_pk_fixed, &
                             sm_pk_Whirl, sm_pk_Shaker
  use Driver_interface, only: Driver_abortFlash
  implicit none

  ! IO
  integer, intent(in) :: ibd
  real, intent(in)    :: time

  ! Internal variables
  type(sm_structure),  pointer :: body
  type(sm_pk_dataset), pointer :: pk
  integer :: nfix, irest, A, i, idx
  real, allocatable, dimension(:,:) :: Xi, v, vd, vdd

  ! Apply kinematics 1 only to entire surface
  irest = 1

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Get the library of kinematics for surface 1
  pk => sm_pk_info( body%restraints_surf(irest)%kinematics_idx )

  ! set the number of nodes (not dofs) that are being fixed
  nfix = body%nnp

  ! Get the mesh (x,y,z) Xi
  allocate( Xi(NDIM, nfix), v(NDIM, nfix), vd(NDIM,nfix), vdd(NDIM,nfix) )
  do A = 1,nfix
     Xi(1,A) = body%x(A)
     Xi(2,A) = body%y(A)
#if NDIM == MDIM
     Xi(3,A) = body%z(A)
#endif
  end do

  !***************************************
  !* Compute kinematics
  !*
  select case( pk%flag  )

  case( SM_PK_FIXED )
     call sm_pk_fixed(time, nfix, Xi, v, vd, vdd, pk%params )

  case( SM_PK_BERMANWANG )
     call sm_pk_Flapping_BermanWang(time, nfix, &
          Xi,v,vd,vdd, pk%params )

  case( SM_PK_WHIRL )
     call sm_pk_Whirl(time, nfix, Xi, v, vd, vdd, pk%params )

  case( SM_PK_SHAKER )
     call sm_pk_Shaker(time, nfix, Xi, v, vd, vdd, pk%params )

  case default
     call Driver_abortFlash('kinematics type is unknown.')

  end select

  ! Apply kinematics to dofs
  do A = 1,nfix
     do i = 1,body%max_dofs_per_node
        idx           = body%ID(i,A)
        body%qn(idx)  = v(i,A)
        body%qdn(idx) = vd(i,A)
        body%qddn(idx)= vdd(i,A)
     end do
  end do

  ! deallocate
  deallocate(Xi, v, vd, vdd)

end subroutine sm_pk_apply_entireBody
