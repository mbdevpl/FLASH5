!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_assembleFluidForce_rigid
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

!#include "Flash.h"
!#include "SolidMechanics.h"
!#include "constants.h"
!#include "Flash_mpi.h"
!#define WRITEFORCE 1


subroutine sm_surf_assembleFluidForce_rigid(ibd)

  use SolidMechanics_data, only: sm_meshMe,sm_NumBodies,sm_bodyInfo, sm_structure
  use Driver_interface, only: Driver_abortFlash, Driver_getDT
  use sm_surf_interface, only : sm_surf_assembleFluidForce_toPoint
  use Driver_data, only : dr_nBegin,dr_nstep,dr_simtime,dr_dt
  use sm_Misc_interface, only: sm_crossProd
  use gr_sbData, ONLY : gr_sbBodyInfo
  use sm_integinterface, ONLY : sm_pc_isPredictor 

  implicit none
#include "SolidMechanics.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: ibd
  
  ! Internal variables
  type(sm_structure),  pointer :: body
  real, dimension(NDIM) :: point, force_pres, force_visc, moment_pres, moment_visc
  real, allocatable, dimension(:) :: qn_master,qdn_master,qddn_master
  real :: xm(MDIM)
  real, dimension(NDIM) :: rhog,xrhog,rhou,xrhou

  integer :: maxdofs
  integer :: i, idx, i_dim, imaster

  character(len=6) :: str_ibd

  logical, save :: first_call      = .true.
  logical, save :: first_call_intu = .true.

  real, allocatable, dimension(:) :: sendbuf
  integer :: ierr
  
  integer :: p
  real, dimension(NPART_PROPS) :: particle
  real, dimension(NDIM) :: f_surf, pos, rpcm
  real, dimension(NDIM) :: force_vol, moment_vol, force_surf, moment_surf
  real :: vlag

  real :: a1,a0,am1
  integer :: flag

  ! Timestep, as all Solidmechanics -> use AB2 in fluid for FSI.
  real :: ibd_dt

  ! Get timestep
  call Driver_getDT(ibd_dt)

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

  force_vol   = 0.
  moment_vol  = 0.
  force_surf  = 0.
  moment_surf = 0.
  rpcm        = 0.

  if (first_call) then
    first_call = .false.
    ! Either Hs_pres, Hs_visc are zero or have been read from checkpoint file. 
    ! Defer.

  else

    ! Actual state in qn, qdn, qddn for master Node:
    maxdofs   = body%max_dofs_per_node
    allocate(qn_master(maxdofs),qdn_master(maxdofs),qddn_master(maxdofs))
    qn_master(1:maxdofs)  = body%  qn(body%ID(1:maxdofs,imaster))
    qdn_master(1:maxdofs) = body% qdn(body%ID(1:maxdofs,imaster))
    qddn_master(1:maxdofs)= body%qddn(body%ID(1:maxdofs,imaster))

    ! Body frame origin actual location:
    point(1:NDIM) = xm(1:NDIM) + qn_master(1:NDIM)

#ifdef FORCES_FROM_VOL_INTEGRALS /* Integrals of Volume + Surface (Uhlmann JCP 2005) */

    write(*,*) 'FORCES_FROM_VOL_INTEGRALS'

    if (first_call_intu) then
      ! Backward Euler:
      a1 = 1.
      a0 =-1.    
      am1= 0.

      if((dr_nstep-dr_nBegin) .gt. CONSTANT_TWO) first_call_intu = .false.

    else
      ! BDF2:
      a1 = 1.5
      a0 =-2.0
      am1= 0.5
    endif 

    ! Hs_pres is the force container:
    body%Hs_pres = 0.
    body%Hs_visc = 0.


    ! Forces due to volume integrals:
    force_vol(1:MDIM) = ( a1*body%int_rhou(1:MDIM,CONSTANT_ONE) +      & ! acceleration
                          a0*body%int_rhou(1:MDIM,CONSTANT_ZERO)+      & ! Use BDF2    
                          am1*body%int_rhou(1:MDIM,-1))/ibd_dt -       & ! ...
                          body%int_rhog(1:MDIM)                          ! gravity
    
    ! Moments due to volume integrals:
    moment_vol(1:MDIM)= ( a1*body%int_xrhou(1:MDIM,CONSTANT_ONE) +     & ! acceleration
                          a0*body%int_xrhou(1:MDIM,CONSTANT_ZERO)+     & ! Use BDF2
                          am1*body%int_xrhou(1:MDIM,-1))/ibd_dt -      & ! ...
                          body%int_xrhog(1:MDIM)                         ! gravity

    call sm_pc_isPredictor(ibd,flag)

    if (flag .eq. SM_TRUE) then 
       ! Here reassign volume integrals:
       sm_BodyInfo(ibd)%int_rhou(:,-1) = sm_BodyInfo(ibd)%int_rhou(:,CONSTANT_ZERO)
       sm_BodyInfo(ibd)%int_rhou(:,CONSTANT_ZERO) = sm_BodyInfo(ibd)%int_rhou(:,CONSTANT_ONE)
       sm_BodyInfo(ibd)%int_xrhou(:,-1)=sm_BodyInfo(ibd)%int_xrhou(:,CONSTANT_ZERO)
       sm_BodyInfo(ibd)%int_xrhou(:,CONSTANT_ZERO)=sm_BodyInfo(ibd)%int_xrhou(:,CONSTANT_ONE)
    endif

    ! Forces from Lagrangian Particles (Surface integrals):
    do p = 1 , gr_sbBodyInfo(ibd)%totalPart

      ! get patch info
      particle = gr_sbBodyInfo(ibd)% particles(1:NPART_PROPS,p)

      vlag = particle(AREA_PART_PROP)*particle(HL_PART_PROP)

      ! Forces:
#if NDIM == 2
      f_surf = (/ particle(FUL_PART_PROP), particle(FVL_PART_PROP) /)
#elif NDIM == 3
      f_surf = (/ particle(FUL_PART_PROP), particle(FVL_PART_PROP),particle(FWL_PART_PROP) /)
#endif

      force_surf(1:NDIM) = force_surf(1:NDIM) + f_surf(1:NDIM)*vlag     

      ! Moments:
      ! position of center of patch
#if NDIM == 2
      pos = (/ particle(POSX_PART_PROP), particle(POSY_PART_PROP) /)
#elif NDIM == 3
      pos = (/ particle(POSX_PART_PROP), particle(POSY_PART_PROP), particle(POSZ_PART_PROP) /)
#endif

     rpcm(1:NDIM) = pos(1:NDIM)-point(1:NDIM)
     ! compute total moment
#if NDIM == 2
      moment_surf(CONSTANT_ONE) = moment_surf(CONSTANT_ONE) + &
                                  sm_crossProd(rpcm,f_surf*vlag)
#elif NDIM == 3
      moment_surf = moment_surf + sm_crossProd(rpcm,f_surf*vlag)
#endif

    enddo

    ! Add to whatever is in Hs: We use Hs_pres as the force container
    ! Forces:
    i_dim = 0
    do i = body%ix,body%ex
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     body%Hs_pres(idx) = body%Hs_pres(idx) + force_vol(i_dim) - force_surf(i_dim)
    end do

    ! Moments: 
    ! Remember, in 2D moment(CONSTANT_ONE) is the moment for theta variable 
    ! (rotation around axis normal to x-y plane): 
    i_dim = 0
    do i = body%iw,body%ew
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     body%Hs_pres(idx) = body%Hs_pres(idx) + moment_vol(i_dim) - moment_surf(i_dim)
    end do

    ! No Metabodies for now.
    if (body%Metabody .gt. CONSTANT_ZERO) &
    call Driver_abortFlash("No Metabodies for FORCES_FROM_VOL_INTEGRALS. Need to implement.")

#else

    ! Hs_pres and Hs_visc to zero:
    body%Hs_pres = 0.
    body%Hs_visc = 0.

    ! Now Call the integration Routine:
    call sm_surf_assembleFluidForce_toPoint(ibd, point, force_pres, force_visc, moment_pres, moment_visc)

    ! Add to whatever is in Hs:
    ! Forces:
    i_dim = 0
    do i = body%ix,body%ex
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     body%Hs_pres(idx) = body%Hs_pres(idx) + force_pres(i_dim) 
     body%Hs_visc(idx) = body%Hs_visc(idx) + force_visc(i_dim)
    end do

    ! Moments: 
    ! Remember, in 2D moment(CONSTANT_ONE) is the moment for theta variable 
    ! (rotation around axis normal to x-y plane): 
    i_dim = 0
    do i = body%iw,body%ew
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     body%Hs_pres(idx) = body%Hs_pres(idx) + moment_pres(i_dim) 
     body%Hs_visc(idx) = body%Hs_visc(idx) + moment_visc(i_dim)
    end do

    ! Here if this body is part of a metabody we will perform an all reduce sum for all 
    ! Pressure and viscous forces. 
    ! All processors managing bodies part of a metabody will run for their particular ibd
    ! through the routines sm_IntegX_advance, sm_assemble_ExtForce, sm_assembleFluidForces_rigid
    ! in this case. At this point they should have their contributions to body%Hs_pres, body%Hs_visc
    ! computed. An all reduce on the group is performed to gather the total corresponding
    ! fluid forces for the set:
    if (body%Metabody .gt. CONSTANT_ZERO) then
      allocate(sendbuf(Body%ndofs))

      ! Pressure forces:
      sendbuf(:) = body%Hs_pres(:)
      call MPI_allReduce(sendbuf, body%Hs_pres, Body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

      ! Viscous forces:
      sendbuf(:) = body%Hs_visc(:)
      call MPI_allReduce(sendbuf, body%Hs_visc, Body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

      deallocate(sendbuf)
    endif

#endif

  !write(*,*) 'AngVel=',qdn_master(7:maxdofs)
  deallocate(qn_master,qdn_master,qddn_master)


  endif


  !if (sm_meshMe .eq. MASTER_PE) then
  !write(*,*) ' '
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_pres(body%ID(body%ix:body%ex,imaster))
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_visc(body%ID(body%ix:body%ex,imaster))
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_pres(body%ID(body%iw:body%ew,imaster))
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_visc(body%ID(body%iw:body%ew,imaster))
  !endif

  ! Add to whatever is in Hs:
  ! Forces:
  do i = body%ix,body%ex
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Hs(idx)  = body%Hs(idx) + body%Hs_pres(idx) + body%Hs_visc(idx)
     endif
  end do

  ! Moments: 
  do i = body%iw,body%ew
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Hs(idx)  = body%Hs(idx) + body%Hs_pres(idx) + body%Hs_visc(idx)
     endif
  end do

#ifdef DEBUG_SOLID
  write(*,*) 'Body, Proc =',ibd,sm_meshMe
!#ifdef FORCES_FROM_VOL_INTEGRALS /* Integrals of Volume + Surface (Uhlmann JCP2005) */
  write(*,*) 'FVol=',force_vol
  write(*,*) 'FLag=',force_surf
  write(*,*) 'FTot=',force_vol-force_surf
  write(*,*) 'MUVol=',moment_vol
  write(*,*) 'MUSrf=',-moment_surf
  write(*,*) 'MUhL=',moment_vol-moment_surf
!#else
  ! Now Call the integration Routine:
  call sm_surf_assembleFluidForce_toPoint(ibd,point,force_pres,force_visc,moment_pres,moment_visc)
  write(*,*) 'FPres=',force_pres
  write(*,*) 'FVisc=',force_visc
  write(*,*) 'FTot=' ,force_pres  +  force_visc
  write(*,*) 'Momt =',moment_pres + moment_visc
!#endif
#endif



  return

end subroutine sm_surf_assembleFluidForce_rigid
