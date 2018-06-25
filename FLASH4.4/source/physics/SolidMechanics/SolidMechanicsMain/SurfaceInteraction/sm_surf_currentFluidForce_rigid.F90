!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_currentFluidForce_rigid
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!! 
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"
#define WRITEFORCE 1


subroutine sm_surf_currentFluidForce_rigid(ibd, ndofs, Hs_pres, Hs_visc)
  use SolidMechanics_data, only: sm_meshMe,sm_NumBodies,sm_bodyInfo, sm_structure
  use Driver_interface, only: Driver_abortFlash
  use sm_surf_interface, only : sm_surf_assembleFluidForce_toPoint
  use Driver_data, only : dr_simtime,dr_dt

  implicit none
  integer, intent(in)  :: ibd, ndofs
  real,    intent(out) :: Hs_pres(ndofs), Hs_visc(ndofs)
  
  ! Internal variables
  type(sm_structure),  pointer :: body
  real, dimension(NDIM) :: point, force_pres, force_visc, moment_pres, moment_visc
  real, allocatable, dimension(:) :: qn_master,qdn_master,qddn_master
  real :: xm(MDIM)

  integer :: maxdofs
  integer :: i, idx, i_dim, imaster

  character(len=6) :: str_ibd

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Extract info from Master:
  ! Reference x,y,z:
  imaster = body % borigin_node
  xm = 0.
  xm(1) = body % x(imaster)
  xm(2) = body % y(imaster)
#if NDIM == MDIM
  xm(3) = body % z(imaster)
#endif          

  force_pres  = 0.
  force_visc  = 0.
  moment_pres = 0.
  moment_visc = 0.

  ! Hs_pres and Hs_visc to zero:
  Hs_pres = 0.
  Hs_visc = 0.

  ! Actual state in qn, qdn, qddn for master Node:
  maxdofs   = body%max_dofs_per_node 
  allocate(qn_master(maxdofs),qdn_master(maxdofs),qddn_master(maxdofs))
  qn_master(1:maxdofs)  = body%  qn(body%ID(1:maxdofs,imaster)) 
  qdn_master(1:maxdofs) = body% qdn(body%ID(1:maxdofs,imaster))
  qddn_master(1:maxdofs)= body%qddn(body%ID(1:maxdofs,imaster))
 
  ! Body frame origin actual location:
  point(1:NDIM) = xm(1:NDIM) + qn_master(1:NDIM)

  ! Now Call the integration Routine:
  call sm_surf_assembleFluidForce_toPoint(ibd, point, force_pres, force_visc, moment_pres, moment_visc)

  ! Add to whatever is in Hs:
  ! Forces:
  i_dim = 0
  do i = body%ix,body%ex
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     Hs_pres(idx) = Hs_pres(idx) + force_pres(i_dim) 
     Hs_visc(idx) = Hs_visc(idx) + force_visc(i_dim)
  end do

  ! Moments: 
  i_dim = 0
  do i = body%iw,body%ew
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     Hs_pres(idx) = Hs_pres(idx) + moment_pres(i_dim) 
     Hs_visc(idx) = Hs_visc(idx) + moment_visc(i_dim)
  end do

#ifdef DEBUG_SOLID
  write(*,*) 'Body, Proc =',ibd,sm_meshMe
  write(*,*) 'FPres=',force_pres
  write(*,*) 'FVisc=',force_visc
  write(*,*) 'Momt =',moment_pres + moment_visc
#endif

  return

end subroutine sm_surf_currentFluidForce_rigid
