!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/forceInBody_analytical/ib_forceInsideBody
!!
!!
!! NAME
!!
!!  ib_forceInsideBody
!!
!!
!! SYNOPSIS
!!
!!  ib_forceInsideBody()
!!
!!
!! DESCRIPTION
!!
!! Routine that forces eulerian velocities inside bodies. Stub.
!! For the moment only analytical shapes are used.
!!
!!***

subroutine ib_forceInsideBody()

  use SolidMechanics_data, only : sm_MeshMe,sm_NumProcs, sm_NumBodies, sm_BodyInfo, sm_meshComm

  use Grid_interface, ONLY : Grid_getListOfBlocks,   &
                             Grid_getDeltas,         &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, only : Driver_abortFlash

  use ImBound_data, only : ib_dt, ib_BlockMarker_flag

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "SolidMechanics.h"


  ! Local Vars
  integer, allocatable, dimension(:) :: flag_forceinside
  integer, allocatable, dimension(:) :: annbody_type,annbody_nparam,aux_v
  real, allocatable, dimension(:,:)  :: annbody_param, r_bod, rd_bod, NwB_Nbod, aux_p,aux_r
  real, allocatable, dimension(:,:,:) :: TNB_bod, aux_t

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  real :: bsize(MDIM),coord(MDIM),del(MDIM)
  integer :: blockCount
  integer, dimension(MAXBLOCKS) :: blockList 
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(2,MDIM) :: boundBox

  real :: xedge(GRID_IHI_GC+1),xcell(GRID_IHI_GC)
  real :: yedge(GRID_JHI_GC+1),ycell(GRID_JHI_GC)
  real :: zedge(GRID_KHI_GC+1),zcell(GRID_KHI_GC)
  
  real :: TBN_xedge(MDIM,GRID_IHI_GC+1),TBN_xcell(MDIM,GRID_IHI_GC)
  real :: TBN_yedge(MDIM,GRID_JHI_GC+1),TBN_ycell(MDIM,GRID_JHI_GC)
  real :: TBN_zedge(MDIM,GRID_KHI_GC+1),TBN_zcell(MDIM,GRID_KHI_GC)

  integer :: i,j,k,ibd,lb,blockID,imaster
  integer :: xpt,ypt,zpt
  real, parameter :: eps = 1.e-12
  real, dimension(MDIM,LOW:HIGH) :: body_box,lb_box
  integer :: pt(MDIM),auxpt(MDIM)
  real :: vec(MDIM),vecnorm
  real :: xvel,yvel,zvel

  real :: TBN(MDIM,MDIM),TNB(MDIM,MDIM)
  real :: locp1(MDIM,1),locp2(MDIM,1),glbp1(MDIM,1),glbp2(MDIM,1)
  real :: radius, cyllen, cyllen2
  integer :: cyldir
  logical :: inlenflg
  integer :: axis_1,axis_2,axis_3

!  integer, save :: uvel_notforced(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS)
!  integer, save :: vvel_notforced(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,MAXBLOCKS)
!  integer, save :: wvel_notforced(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1,MAXBLOCKS)

!  real, save :: uvel_val(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS)
!  real, save :: vvel_val(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,MAXBLOCKS)
!  real, save :: wvel_val(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1,MAXBLOCKS)
  
  integer :: ierr

  logical, save :: firstcall =.true.

  !return

  ! Allocations:
  allocate( flag_forceinside(sm_NumBodies), &
            annbody_type(sm_NumBodies)    , &
            annbody_nparam(sm_NumBodies), aux_v(sm_NumBodies) )
  allocate( annbody_param(RB_MAXANNPARAM,sm_NumBodies), aux_p(RB_MAXANNPARAM,sm_NumBodies) )
  allocate( r_bod(MDIM,sm_NumBodies) , rd_bod(MDIM,sm_NumBodies) , &
            NwB_Nbod(MDIM,sm_NumBodies), aux_r(MDIM,sm_NumBodies) )
  allocate( TNB_bod(MDIM,MDIM,sm_NumBodies), aux_t(MDIM,MDIM,sm_NumBodies) )


!  if (firstcall) then
!    uvel_notforced = 1
!    vvel_notforced = 1
!    wvel_notforced = 1
!    uvel_val = 0.
!    vvel_val = 0.
!    wvel_val = 0.
!  endif

  flag_forceinside = SM_FALSE
  annbody_type     = 0
  annbody_nparam   = 0
  annbody_param    = 0.
  r_bod = 0.; rd_bod = 0.; NwB_Nbod = 0.;
  TNB_bod = 0.

  ! First populate body info
  do ibd = 1,sm_NumBodies
     if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

        ! Parameters:
        flag_forceinside(ibd) = sm_BodyInfo(ibd) % flag_forceinside
        annbody_type(ibd)     = sm_BodyInfo(ibd) % annbody_type
        annbody_nparam(ibd)   = sm_BodyInfo(ibd) % annbody_nparam
        if( annbody_nparam(ibd) .gt. 0 ) &
        annbody_param(1:annbody_nparam(ibd),ibd) = sm_BodyInfo(ibd) % annbody_param(1:annbody_nparam(ibd))

        ! Absolute positions, velocities and angular velocity:
        imaster = sm_BodyInfo(ibd) % borigin_node

        r_bod(IAXIS,ibd)=sm_BodyInfo(ibd)%x(imaster) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(IAXIS,imaster))
        !r_bod(JAXIS,ibd)=sm_BodyInfo(ibd)%x(imaster) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(JAXIS,imaster)) ! the %x should be %y, Shizhao
        r_bod(JAXIS,ibd)=sm_BodyInfo(ibd)%y(imaster) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(JAXIS,imaster))  ! Shizhao
#if NDIM == MDIM
        !r_bod(KAXIS,ibd)=sm_BodyInfo(ibd)%x(imaster) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(KAXIS,imaster)) ! The %x should be %z, Shizhao
        r_bod(KAXIS,ibd)=sm_BodyInfo(ibd)%z(imaster) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(KAXIS,imaster))  ! Shizhao
#else
        r_bod(KAXIS,ibd)=0.0
#endif
        rd_bod(IAXIS,ibd)= sm_BodyInfo(ibd)%qdn(sm_BodyInfo(ibd)%ID(IAXIS,imaster))
        rd_bod(JAXIS,ibd)= sm_BodyInfo(ibd)%qdn(sm_BodyInfo(ibd)%ID(JAXIS,imaster))
#if NDIM == MDIM 
        rd_bod(KAXIS,ibd)= sm_BodyInfo(ibd)%qdn(sm_BodyInfo(ibd)%ID(KAXIS,imaster))
#else
        rd_bod(KAXIS,ibd)=0.0
#endif

        NwB_Nbod(1:MDIM,ibd) = sm_BodyInfo(ibd)%RBMAT%NwB_N(1:MDIM,1)
        

        ! Actual transformation matrix:
        TNB_bod(1:MDIM,1:MDIM,ibd) =  sm_BodyInfo(ibd)%RBMAT%TNB(1:MDIM,1:MDIM)

     endif
  enddo

!  call Timers_start("forceinside_allreduce")
  ! Now All reduce.... ughh expensive.
  ! Parameters
  aux_v = flag_forceinside
  call mpi_allreduce ( aux_v, flag_forceinside, sm_NumBodies, FLASH_INTEGER, &
                       MPI_SUM, sm_meshComm, ierr )
  aux_v = annbody_type
  call mpi_allreduce ( aux_v, annbody_type, sm_NumBodies, FLASH_INTEGER, &
                       MPI_SUM, sm_meshComm, ierr )
  aux_v = annbody_nparam
  call mpi_allreduce ( aux_v, annbody_nparam, sm_NumBodies, FLASH_INTEGER, &
                       MPI_SUM, sm_meshComm, ierr )

  ! Parameter list:
  aux_p = annbody_param 
  call mpi_allreduce ( aux_p, annbody_param,RB_MAXANNPARAM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr )                   

  ! Positions Velocities and ang vels:
  aux_r = r_bod
  call mpi_allreduce ( aux_r, r_bod,MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr ) 

  aux_r = rd_bod
  call mpi_allreduce ( aux_r, rd_bod,MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr ) 

  aux_r = NwB_Nbod
  call mpi_allreduce ( aux_r, NwB_Nbod,MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr ) 

  ! TNB:
  aux_t = TNB_bod
  call mpi_allreduce ( aux_t, TNB_bod,MDIM*MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr ) 
!  call Timers_stop("forceinside_allreduce")

  ! Then, Block by block run over bodies to force internal eulerian points:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  pt = 0;
  ! Now do loop over Grid Blocks: 
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get Blocks center coords:
     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get Blocks bsize:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Block box including guardcells:
     lb_box(IAXIS,LOW)  = coord(IAXIS) - bsize(IAXIS)/2.0 - real(NGUARD)*del(IAXIS)
     lb_box(IAXIS,HIGH) = coord(IAXIS) + bsize(IAXIS)/2.0 + real(NGUARD)*del(IAXIS)
     lb_box(JAXIS,LOW)  = coord(JAXIS) - bsize(JAXIS)/2.0 - real(NGUARD)*del(JAXIS)
     lb_box(JAXIS,HIGH) = coord(JAXIS) + bsize(JAXIS)/2.0 + real(NGUARD)*del(JAXIS)        
#if NDIM == MDIM
     lb_box(KAXIS,LOW)  = coord(KAXIS) - bsize(KAXIS)/2.0 - real(NGUARD)*del(KAXIS)
     lb_box(KAXIS,HIGH) = coord(KAXIS) + bsize(KAXIS)/2.0 + real(NGUARD)*del(KAXIS)
#else
     lb_box(KAXIS,LOW)  = 0.0
     lb_box(KAXIS,HIGH) = 0.0
#endif
     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        
     ! Compute Grid line locations: 
#if NDIM == MDIM
     ! Z:
     do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        zedge(k) = coord(KAXIS) - bsize(KAXIS)/2.0 + real(k - NGUARD - 1)*del(KAXIS)
        zcell(k) = zedge(k) + 0.5*del(KAXIS)     
     enddo 
     zedge(blkLimitsGC(HIGH,KAXIS)+1)=coord(KAXIS)-bsize(KAXIS)/2.0 + &
                                      real(blkLimitsGC(HIGH,KAXIS)-NGUARD)*del(KAXIS)
#else
     zedge(:) = 0.0
     zcell(:) = 0.0
#endif

     ! Y:
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        yedge(j) = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j - NGUARD - 1)*del(JAXIS)
        ycell(j) = yedge(j) + 0.5*del(JAXIS)
     enddo
     yedge(blkLimitsGC(HIGH,JAXIS)+1)=coord(JAXIS)-bsize(JAXIS)/2.0 + &
                                      real(blkLimitsGC(HIGH,JAXIS)-NGUARD)*del(JAXIS)

     ! X:
     do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
        xedge(i) = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i - NGUARD - 1)*del(IAXIS)
        xcell(i) = xedge(i) + 0.5*del(IAXIS)
     enddo
     xedge(blkLimitsGC(HIGH,IAXIS)+1)=coord(IAXIS)-bsize(IAXIS)/2.0+ &
                                      real(blkLimitsGC(HIGH,IAXIS)-NGUARD)*del(IAXIS)

     ! Point to blocks face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! Now loop over bodies:
     do ibd = 1,sm_NumBodies

        if( flag_forceinside(ibd) .eq. SM_FALSE ) cycle

        ! TBN matrix: TBN = transpose( TNB_bod(1:MDIM,1:MDIM,ibd) )
        do j=1,MDIM
           do i=1,MDIM
              TBN(i,j) = TNB_bod(j,i,ibd)
           enddo
        enddo

        ! First individualize bodies bounding box:
        select case ( annbody_type(ibd) )

        case ( RB_ANNSPHERE )

#if NDIM == MDIM              
           ! Sphere:
           radius = annbody_param(1,ibd)/2.         
           body_box(IAXIS,LOW)  = r_bod(IAXIS,ibd) - radius
           body_box(IAXIS,HIGH) = r_bod(IAXIS,ibd) + radius
           body_box(JAXIS,LOW)  = r_bod(JAXIS,ibd) - radius
           body_box(JAXIS,HIGH) = r_bod(JAXIS,ibd) + radius
           body_box(KAXIS,LOW)  = r_bod(KAXIS,ibd) - radius
           body_box(KAXIS,HIGH) = r_bod(KAXIS,ibd) + radius
           
#else
           call Driver_abortFlash("Cannot Force sphere in 2D. Use RB_ANNDISC")

#endif           
        case ( RB_ANNDISC )
           
           cyldir = int(annbody_param(1,ibd))
           radius = annbody_param(2,ibd)/2.
           cyllen2= annbody_param(3,ibd)/2.

           TNB(1:MDIM,1:MDIM) =  TNB_bod(1:MDIM,1:MDIM,ibd)

           ! Local xylinder axis extrema locations:
           locp1(1:MDIM,1) = 0.; locp2(1:MDIM,1) = 0.;
           select case (cyldir)
           case(IAXIS) ! Cylinder direction in body axis xB  
 
              ! Local position of lower cylinder center:  (p1) x----cen----x (p2)
              locp1(IAXIS,1) = -cyllen2

              ! Local position of higher cylinder center:
              locp2(IAXIS,1) =  cyllen2

           case(JAXIS) ! Cylinder direction in body axis yB  

              ! Local position of lower cylinder center:  (p1) x----cen----x (p2)
              locp1(JAXIS,1) = -cyllen2

              ! Local position of higher cylinder center:
              locp2(JAXIS,1) =  cyllen2

           case(KAXIS)

              ! Local position of lower cylinder center:  (p1) x----cen----x (p2)
              locp1(KAXIS,1) = -cyllen2

              ! Local position of higher cylinder center:
              locp2(KAXIS,1) =  cyllen2

           end select

           ! Global positions of extreme points p1,p2:  
           glbp1 = matmul( TNB , locp1 )  
           glbp2 = matmul( TNB , locp2 )              

           ! Approx body_box: (It adds radius in each dir, which is an overkill)
#if NDIM == MDIM
           body_box(IAXIS,LOW)  = r_bod(IAXIS,ibd) + &
                                  min(glbp1(IAXIS,1),glbp2(IAXIS,1)) - radius
           body_box(IAXIS,HIGH) = r_bod(IAXIS,ibd) + &
                                  max(glbp1(IAXIS,1),glbp2(IAXIS,1)) + radius
           body_box(JAXIS,LOW)  = r_bod(JAXIS,ibd) + &
                                  min(glbp1(JAXIS,1),glbp2(JAXIS,1)) - radius
           body_box(JAXIS,HIGH) = r_bod(JAXIS,ibd) + &
                                  max(glbp1(JAXIS,1),glbp2(JAXIS,1)) + radius
           body_box(KAXIS,LOW)  = r_bod(KAXIS,ibd) + &
                                  min(glbp1(KAXIS,1),glbp2(KAXIS,1)) - radius
           body_box(KAXIS,HIGH) = r_bod(KAXIS,ibd) + &
                                  max(glbp1(KAXIS,1),glbp2(KAXIS,1)) + radius
#else
           body_box(IAXIS,LOW)  = r_bod(IAXIS,ibd) - radius
           body_box(IAXIS,HIGH) = r_bod(IAXIS,ibd) + radius
           body_box(JAXIS,LOW)  = r_bod(JAXIS,ibd) - radius
           body_box(JAXIS,HIGH) = r_bod(JAXIS,ibd) + radius           
           body_box(KAXIS,LOW)  = 0.0
           body_box(KAXIS,HIGH) = 0.0
#endif

!           call Driver_abortFlash('ib_forceInsideBody: Disc analytical geometry operations not coded yet.')
           
        case ( RB_ANNRBC  )
           
           call Driver_abortFlash('ib_forceInsideBody: RBC analytical geometry operations not coded yet.')
           
        case default
           
           call Driver_abortFlash('ib_forceInsideBody: Analytical geometry type not known.')
           
        end select

        ! Test to see if there is no overlap:
        if ((body_box(IAXIS,HIGH).lt. lb_box(IAXIS,LOW)) .or. &
            (body_box(IAXIS,LOW) .gt. lb_box(IAXIS,HIGH))) cycle
        if ((body_box(JAXIS,HIGH).lt. lb_box(JAXIS,LOW)) .or. &
            (body_box(JAXIS,LOW) .gt. lb_box(JAXIS,HIGH))) cycle
#if NDIM == MDIM
        if ((body_box(KAXIS,HIGH).lt. lb_box(KAXIS,LOW)) .or. &
            (body_box(KAXIS,LOW) .gt. lb_box(KAXIS,HIGH))) cycle
#endif

        ! Transform Grid line coordinates to Body Frame: 
#if NDIM == MDIM
        ! Z:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           TBN_zedge(1:MDIM,k) = (zedge(k)-r_bod(KAXIS,ibd)) * TBN(1:MDIM,KAXIS)
           TBN_zcell(1:MDIM,k) = (zcell(k)-r_bod(KAXIS,ibd)) * TBN(1:MDIM,KAXIS)    
        enddo
        TBN_zedge(1:MDIM,blkLimitsGC(HIGH,KAXIS)+1) = &
                 (zedge(blkLimitsGC(HIGH,KAXIS)+1)-r_bod(KAXIS,ibd)) * TBN(1:MDIM,KAXIS)
#else
        TBN_zedge(:,:) = 0.0
        TBN_zcell(:,:) = 0.0
#endif

        ! Y:
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           TBN_yedge(1:MDIM,j) = (yedge(j)-r_bod(JAXIS,ibd)) * TBN(1:MDIM,JAXIS)
           TBN_ycell(1:MDIM,j) = (ycell(j)-r_bod(JAXIS,ibd)) * TBN(1:MDIM,JAXIS)
        enddo
        TBN_yedge(1:MDIM,blkLimitsGC(HIGH,JAXIS)+1) = &
                 (yedge(blkLimitsGC(HIGH,JAXIS)+1)-r_bod(JAXIS,ibd)) * TBN(1:MDIM,JAXIS)

        ! X:
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           TBN_xedge(1:MDIM,i) = (xedge(i)-r_bod(IAXIS,ibd)) * TBN(1:MDIM,IAXIS)
           TBN_xcell(1:MDIM,i) = (xcell(i)-r_bod(IAXIS,ibd)) * TBN(1:MDIM,IAXIS)
        enddo
        TBN_xedge(1:MDIM,blkLimitsGC(HIGH,IAXIS)+1) = &
                 (xedge(blkLimitsGC(HIGH,IAXIS)+1)-r_bod(IAXIS,ibd)) * TBN(1:MDIM,IAXIS)


        ! Loop over points, test if they belong to body and force:
        select case ( annbody_type(ibd) ) 
        case ( RB_ANNSPHERE )

#if NDIM == MDIM
           radius = annbody_param(1,ibd)/2.  

           ! Z velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                    if (abs(facezData(FORC_FACE_VAR,i,j,k)) .lt. eps) then
                    !if ((wvel_notforced(i,j,k,blockID) .eq. 1) .and. (abs(facezData(FORC_FACE_VAR,i,j,k)) .lt. eps)) then
                       ! Distance from body center of mass to Z veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zedge(1:MDIM,k)
                       vecnorm = sqrt(vec(IAXIS)**2. + vec(JAXIS)**2. + vec(KAXIS)**2.)

                       ! Inside sphere radius?
                       if ( vecnorm .lt. radius ) then
                          ! vz = vzo + (wx*ypt - wy*xpt)
                          zvel = rd_bod(KAXIS,ibd) +  NwB_Nbod(IAXIS,ibd)*(ycell(j)-r_bod(JAXIS,ibd)) - &
                                                      NwB_Nbod(JAXIS,ibd)*(xcell(i)-r_bod(IAXIS,ibd))
                          facezData(FORC_FACE_VAR,i,j,k) = facezData(FORC_FACE_VAR,i,j,k) + &
                                                           (zvel-facezData(VELC_FACE_VAR,i,j,k))/ib_dt
                          !wvel_notforced(i,j,k,blockID) = 0
                          !wvel_val(i,j,k,blockID) = zvel
                          if (k .ne. (blkLimits(HIGH,KAXIS)+1) ) pt(KAXIS) = pt(KAXIS) + 1
                       endif
                    endif

                 enddo
              enddo
           enddo

           ! Y velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                    if (abs(faceyData(FORC_FACE_VAR,i,j,k)) .lt. eps) then
                    !if ((vvel_notforced(i,j,k,blockID) .eq. 1) .and. (abs(faceyData(FORC_FACE_VAR,i,j,k)) .lt. eps)) then
                       ! Distance from body center of mass to Y veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_yedge(1:MDIM,j) + TBN_zcell(1:MDIM,k)
                       vecnorm = sqrt(vec(IAXIS)**2. + vec(JAXIS)**2. + vec(KAXIS)**2.)

                       ! Inside sphere radius?
                       if ( vecnorm .lt. radius ) then
                          ! vy = vyo + (wz*xpt - wx*zpt)
                          yvel = rd_bod(JAXIS,ibd) +  NwB_Nbod(KAXIS,ibd)*(xcell(i)-r_bod(IAXIS,ibd)) - &
                                                      NwB_Nbod(IAXIS,ibd)*(zcell(k)-r_bod(KAXIS,ibd))
                          faceyData(FORC_FACE_VAR,i,j,k) = faceyData(FORC_FACE_VAR,i,j,k) + &
                                                           (yvel-faceyData(VELC_FACE_VAR,i,j,k))/ib_dt
                          !vvel_notforced(i,j,k,blockID) = 0
                          !vvel_val(i,j,k,blockID) = yvel
                          if (j .ne. (blkLimits(HIGH,JAXIS)+1)) pt(JAXIS) = pt(JAXIS) + 1
                       endif
                    endif

                 enddo
              enddo
           enddo

           ! X velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1

                    if (abs(facexData(FORC_FACE_VAR,i,j,k)) .lt. eps) then
                    !if ((uvel_notforced(i,j,k,blockID) .eq. 1) .and. (abs(facexData(FORC_FACE_VAR,i,j,k)) .lt. eps)) then
                       ! Distance from body center of mass to X veloc point
                       vec(1:MDIM) = TBN_xedge(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zcell(1:MDIM,k)
                       vecnorm = sqrt(vec(IAXIS)**2. + vec(JAXIS)**2. + vec(KAXIS)**2.)

                       ! Inside sphere radius?
                       if ( vecnorm .lt. radius ) then
                          ! vx = vxo + (wy*zpt - wz*ypt)
                          xvel = rd_bod(IAXIS,ibd) +  NwB_Nbod(JAXIS,ibd)*(zcell(k)-r_bod(KAXIS,ibd)) - &
                                                      NwB_Nbod(KAXIS,ibd)*(ycell(j)-r_bod(JAXIS,ibd))
                          facexData(FORC_FACE_VAR,i,j,k) = facexData(FORC_FACE_VAR,i,j,k) + &
                                                           (xvel-facexData(VELC_FACE_VAR,i,j,k))/ib_dt
                          !uvel_notforced(i,j,k,blockID) = 0
                          !uvel_val(i,j,k,blockID) = xvel
                          if (i .ne. (blkLimits(HIGH,IAXIS)+1)) pt(IAXIS) = pt(IAXIS) + 1
                       endif
                    endif

                 enddo
              enddo
           enddo
#else
           call Driver_abortFlash("Cannot Force sphere in 2D. Use RB_ANNDISC")

#endif 
        case ( RB_ANNDISC )

           cyldir = int(annbody_param(1,ibd))
           radius = annbody_param(2,ibd)/2. + eps
           cyllen = annbody_param(3,ibd)
           cyllen2= cyllen/2. + eps

           select case (cyldir)
           case(IAXIS) ! Cylinder direction in body axis xB  

              axis_1 = IAXIS
              axis_2 = JAXIS
              axis_3 = KAXIS

           case(JAXIS) ! Cylinder direction in body axis yB  

              axis_1 = JAXIS
              axis_2 = KAXIS
              axis_3 = IAXIS

           case(KAXIS) ! Cylinder direction in body axis zB

              axis_1 = KAXIS
              axis_2 = IAXIS
              axis_3 = JAXIS

           end select

           call Timers_start("velocity_insideforce") 

#if NDIM == MDIM
           ! Z velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                    if (abs(facezData(FORC_FACE_VAR,i,j,k)) .lt. eps) then
                    !if ((wvel_notforced(i,j,k,blockID) .eq. 1) .and. (abs(facezData(FORC_FACE_VAR,i,j,k)) .lt. eps)) then
                       ! Distance from body center of mass to Z veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zedge(1:MDIM,k)

                       inlenflg= abs(vec(axis_1)) .lt. (cyllen2)
                       vecnorm = sqrt(vec(axis_2)**2. + vec(axis_3)**2.)

                       ! Inside cylinder?
                       if ( inlenflg .and. (vecnorm .lt. (radius)) ) then
                          ! vz = vzo + (wx*ypt - wy*xpt)
                          zvel = rd_bod(KAXIS,ibd) +  NwB_Nbod(IAXIS,ibd)*(ycell(j)-r_bod(JAXIS,ibd)) - &
                                                      NwB_Nbod(JAXIS,ibd)*(xcell(i)-r_bod(IAXIS,ibd))
                          facezData(FORC_FACE_VAR,i,j,k) = facezData(FORC_FACE_VAR,i,j,k) + &
                                                           (zvel-facezData(VELC_FACE_VAR,i,j,k))/ib_dt
                          !wvel_notforced(i,j,k,blockID)  = 0
                          !wvel_val(i,j,k,blockID)  = zvel
                          if (k .lt. (blkLimits(HIGH,KAXIS)+1)) pt(KAXIS) = pt(KAXIS) + 1
                       endif
                    endif

                 enddo
              enddo
           enddo
#endif
           ! Y velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                    if (abs(faceyData(FORC_FACE_VAR,i,j,k)) .lt. eps) then
                    !if ((vvel_notforced(i,j,k,blockID) .eq. 1) .and. (abs(faceyData(FORC_FACE_VAR,i,j,k)) .lt. eps)) then
                       ! Distance from body center of mass to Y veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_yedge(1:MDIM,j) + TBN_zcell(1:MDIM,k)

                       inlenflg= abs(vec(axis_1)) .lt. (cyllen2)
                       vecnorm = sqrt(vec(axis_2)**2. + vec(axis_3)**2.)

                       ! Inside cylinder?
                       if ( inlenflg .and. (vecnorm .lt. (radius)) ) then
#if NDIM == MDIM
                          ! vy = vyo + (wz*xpt - wx*zpt)
                          yvel = rd_bod(JAXIS,ibd) +  NwB_Nbod(KAXIS,ibd)*(xcell(i)-r_bod(IAXIS,ibd)) - &
                                                      NwB_Nbod(IAXIS,ibd)*(zcell(k)-r_bod(KAXIS,ibd))
#else
                          ! vy = vyo + wz*xpt, wz in location 1
                          yvel = rd_bod(JAXIS,ibd) + NwB_Nbod(1,ibd)*(xcell(i)-r_bod(IAXIS,ibd))
#endif
                          faceyData(FORC_FACE_VAR,i,j,k) = faceyData(FORC_FACE_VAR,i,j,k) + &
                                                           (yvel-faceyData(VELC_FACE_VAR,i,j,k))/ib_dt
                          !vvel_notforced(i,j,k,blockID)  = 0
                          !vvel_val(i,j,k,blockID)  = yvel
                          if (j .lt. (blkLimits(HIGH,JAXIS)+1)) pt(JAXIS) = pt(JAXIS) + 1
                       endif
                    endif

                 enddo
              enddo
           enddo
           ! X velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1

                    if (abs(facexData(FORC_FACE_VAR,i,j,k)) .lt. eps) then
                    !if ((uvel_notforced(i,j,k,blockID) .eq. 1) .and. (abs(facexData(FORC_FACE_VAR,i,j,k)) .lt. eps)) then
                       ! Distance from body center of mass to X veloc point
                       vec(1:MDIM) = TBN_xedge(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zcell(1:MDIM,k)

                       inlenflg= abs(vec(axis_1)) .lt. (cyllen2)
                       vecnorm = sqrt(vec(axis_2)**2. + vec(axis_3)**2.)

                       ! Inside cylinder?
                       if ( inlenflg .and. (vecnorm .lt. (radius)) ) then
#if NDIM == MDIM
                          ! vx = vxo + (wy*zpt - wz*ypt)
                          xvel = rd_bod(IAXIS,ibd) +  NwB_Nbod(JAXIS,ibd)*(zcell(k)-r_bod(KAXIS,ibd)) - &
                                                      NwB_Nbod(KAXIS,ibd)*(ycell(j)-r_bod(JAXIS,ibd))
#else
                          ! vx = vxo - wz*ypt, wz in location 1
                          xvel = rd_bod(IAXIS,ibd) - NwB_Nbod(1,ibd)*(ycell(j)-r_bod(JAXIS,ibd))
#endif
                          facexData(FORC_FACE_VAR,i,j,k) = facexData(FORC_FACE_VAR,i,j,k) + &
                                                           (xvel-facexData(VELC_FACE_VAR,i,j,k))/ib_dt
                          !uvel_notforced(i,j,k,blockID)  = 0
                          !uvel_val(i,j,k,blockID)  = xvel
                          if (i .lt. (blkLimits(HIGH,IAXIS)+1)) pt(IAXIS) = pt(IAXIS) + 1
                       endif
                    endif

                 enddo
              enddo
           enddo
           call Timers_stop("velocity_insideforce")

           
        case ( RB_ANNRBC  )
           
           call Driver_abortFlash('ib_forceInsideBody: RBC analytical geometry operations not coded yet.')
           
        case default
           
           call Driver_abortFlash('ib_forceInsideBody: Analytical geometry type not known.')
           
        end select


     enddo

     ! Release Pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo
  
  ! Total Points in u,v,w forced:
  auxpt = pt
  call mpi_allreduce( auxpt, pt, MDIM, FLASH_INTEGER, MPI_SUM, sm_meshComm, ierr )

! Turn off to reduce the time in printing on screens. Shizhao Wang, Jan 07, 2015.
!  if (sm_meshMe .eq. MASTER_PE) write(*,*) 'Total Internal Points forced u,v,w=',pt(IAXIS:KAXIS)  

  ! Deallocate:
  deallocate( flag_forceinside ,  annbody_type , annbody_nparam )
  deallocate( annbody_param )
  deallocate( r_bod , rd_bod , NwB_Nbod )
  deallocate( TNB_bod )
  deallocate( aux_v , aux_p , aux_r , aux_t ) 

  return

end subroutine ib_forceInsideBody
