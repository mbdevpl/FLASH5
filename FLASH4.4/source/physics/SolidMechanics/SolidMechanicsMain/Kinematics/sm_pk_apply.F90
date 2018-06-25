!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_pk_apply
!!
!!
!! NAME
!!
!! 
!!
!!
!! SYNOPSIS
!! 
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
#include "SolidMechanics.h"

subroutine sm_pk_apply(ibd,time)

  use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
  use sm_pk_data, only: sm_pk_dataset, sm_pk_info, sm_pk_timedelay
  use sm_pk_interface, only: sm_pk_Flapping_BermanWang, sm_pk_fixed, sm_pk_Whirl, &
                             sm_pk_fixed_dof,sm_pk_harmonic_dof, &
                             sm_pk_Shaker, sm_pk_constvel_dof
  use Driver_interface, only: Driver_abortFlash
  implicit none

  ! IO
  integer, intent(in) :: ibd
  real, intent(in)    :: time
  real, allocatable, dimension(:) :: vc, vcd, vcdd, paramcoord
  integer :: nfix_coord,maxrestparams,ircoord

  ! Internal variables
  type(sm_structure),  pointer :: body
  type(sm_pk_dataset), pointer :: pk
  integer :: nfix, irest, A, i, j, idx, maxdofs
  real, allocatable, dimension(:,:) :: Xi, v, vd, vdd
  real :: time_mod

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! shift the time per the specified delay
  time_mod = time - sm_pk_timedelay

  ! Set Restrained individual nodes for Body:
  nfix_coord = Body%nrcoords
  if (nfix_coord .gt. 0) then 

     allocate( vc(nfix_coord), vcd(nfix_coord), vcdd(nfix_coord) )
     maxrestparams = body%maxrestparams
     allocate( paramcoord(maxrestparams) )

     do ircoord=1,nfix_coord

        paramcoord(1:maxrestparams) = body%restraints(ircoord)%param(1:maxrestparams)

        select case( body%restraints(ircoord)%restype )
        case(SM_PK_FIXED)

           call sm_pk_fixed_dof(time_mod,maxrestparams,paramcoord,vc(ircoord),vcd(ircoord),vcdd(ircoord))

        case(SM_PK_CONSTVEL)

           
           call sm_pk_constvel_dof(time_mod,maxrestparams,paramcoord,vc(ircoord),vcd(ircoord),vcdd(ircoord))

        case(SM_PK_HARMONIC)

           call sm_pk_harmonic_dof(time_mod,maxrestparams,paramcoord,vc(ircoord),vcd(ircoord),vcdd(ircoord))
            
        ! MORE CASES HERE ---->
 
        case default
           call Driver_abortFlash('Restrained Dofs: kinematics type is unknown.')
               
        end select
     enddo


     ! Apply kinematics to dofs
     do ircoord = 1,nfix_coord
        A = body%restraints(ircoord)%restnode
        i = body%restraints(ircoord)%restdof
        idx           = body%ID(i,A)
        body%qn(idx)  = vc(ircoord)
        body%qdn(idx) = vcd(ircoord)
        body%qddn(idx)= vcdd(ircoord)
     end do
      
     ! deallocate
     deallocate(vc, vcd, vcdd, paramcoord)

  endif


  ! surfaces
  do irest = 1,body%nrsurf

     ! Get the library of kinematics
     pk => sm_pk_info( body%restraints_surf(irest)%kinematics_idx )

     ! set the number of nodes (not dofs) that are being fixed
     nfix = body%restraints_surf(irest)%nfix

     ! Get the mesh (x,y,z) Xi
     maxdofs   = body%max_dofs_per_node
     allocate( Xi(MDIM,nfix), v(maxdofs,nfix), vd(maxdofs,nfix), vdd(maxdofs,nfix) )
     Xi = 0.
     do j = 1,nfix
        A = body%restraints_surf(irest)%node_list(j)
        Xi(1,j) = body%x(A)
        Xi(2,j) = body%y(A)
#if NDIM == MDIM
        Xi(3,j) = body%z(A)
#endif
     end do

     !***************************************
     !* Compute kinematics
     !*
     select case( pk%flag  )

     case( SM_PK_FIXED )
        call sm_pk_fixed(time_mod, nfix, Xi, v, vd, vdd, pk%params )

     case( SM_PK_BERMANWANG )
        call sm_pk_Flapping_BermanWang(time_mod, nfix, &
                                       Xi,v,vd,vdd, pk%params)

     case( SM_PK_WHIRL )
        call sm_pk_Whirl(time_mod, nfix, Xi, v, vd, vdd, pk%params )

     case( SM_PK_SHAKER )
        call sm_pk_Shaker(time_mod, nfix, Xi, v, vd, vdd, pk%params )

     case( SM_PK_MASTERSLAVE )
        ! Postponed.
        goto 163 ! Infamous s.t. no garbage is put in qn, qdn, qddn !!!

     case default
        call Driver_abortFlash('kinematics type is unknown.')

     end select

     ! Apply kinematics to dofs
     if( body%restraints_surf(irest)%restdof == ALL_DOF ) then
        ! All to all the DOFS on the node (i = 1,2,3)
        do j = 1,nfix
           A = body%restraints_surf(irest)%node_list(j)
           do i = 1,body%max_dofs_per_node
              idx           = body%ID(i,A)
              body%qn(idx)  = v(i,j)
              body%qdn(idx) = vd(i,j)
              body%qddn(idx)= vdd(i,j)
           end do
        end do
     else
        ! Apply only specified local coordinate number i
         do j = 1,nfix
           A = body%restraints_surf(irest)%node_list(j)
           i = body%restraints_surf(irest)%restdof
           idx           = body%ID(i,A)
           body%qn(idx)  = v(i,j)
           body%qdn(idx) = vd(i,j)
           body%qddn(idx)= vdd(i,j)
        end do
     end if

     ! deallocate
163  continue
     deallocate(Xi, v, vd, vdd)
     nullify(pk)

  end do

  return 
end subroutine sm_pk_apply
