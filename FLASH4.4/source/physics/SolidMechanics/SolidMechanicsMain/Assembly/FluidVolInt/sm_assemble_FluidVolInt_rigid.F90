! Function to compute fluid integrals on volume occupied by
! rigid bodies.
! Used on Integral form of nodal force computations (Uhlman 2005).


subroutine sm_assemble_FluidVolInt_rigid()

  use SolidMechanics_data, only : sm_MeshMe,sm_NumProcs, sm_NumBodies, sm_BodyInfo, sm_meshComm, &
                                  sm_gravX, sm_gravY, sm_gravZ

  use Grid_interface, ONLY : Grid_getListOfBlocks,   &
                             Grid_getDeltas,         &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, only : Driver_abortFlash

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "SolidMechanics.h"

  integer, allocatable, dimension(:) :: annbody_type,annbody_nparam,aux_v
  real, allocatable, dimension(:,:)  :: annbody_param, r_bod, rd_bod, NwB_Nbod,aux_p,aux_r
  real, allocatable, dimension(:,:,:):: TNB_bod, aux_t

  real, allocatable, dimension(:,:)  :: int_tvol, int_rhog, int_xrhog, int_rhou, int_xrhou 
  real :: rhogx,rhogy,rhogz,rhou,rhov,rhow
  real :: rhoDygx,rhoDzgx,rhoDxgy,rhoDzgy,rhoDxgz,rhoDygz
  real :: rhoDyu,rhoDzu,rhoDxv,rhoDzv,rhoDxw,rhoDyw

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  real :: bsize(MDIM),coord(MDIM),del(MDIM),Vcell,hcell
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
  integer :: pt(MDIM),auxpt(MDIM),pt_damped(MDIM)
  real :: vec(MDIM),vecnorm
  real :: xvel,yvel,zvel

  real :: TBN(MDIM,MDIM),TNB(MDIM,MDIM)
  real :: locp1(MDIM,1),locp2(MDIM,1),glbp1(MDIM,1),glbp2(MDIM,1)
  real :: radius, radius1, radius2, radius3, cyllen, cyllen2
  integer :: cyldir
  logical :: inlenflg
  integer :: axis_1,axis_2,axis_3

  real :: Tvolx, Tvoly, Tvolz, body_vol, vol_fctx, vol_fcty, vol_fctz

  real :: alphaijk, phi(2**NDIM), xc(MDIM,2**NDIM)
  integer :: iph

  integer :: ierr

  real, parameter :: rhoi = 1.

  allocate( annbody_type(sm_NumBodies)    , &
            annbody_nparam(sm_NumBodies), aux_v(sm_NumBodies) )
  allocate( annbody_param(RB_MAXANNPARAM,sm_NumBodies),aux_p(RB_MAXANNPARAM,sm_NumBodies) )
  allocate( r_bod(MDIM,sm_NumBodies) , rd_bod(MDIM,sm_NumBodies) , &
            NwB_Nbod(MDIM,sm_NumBodies), aux_r(MDIM,sm_NumBodies) )
  allocate( TNB_bod(MDIM,MDIM,sm_NumBodies), aux_t(MDIM,MDIM,sm_NumBodies) )

  allocate( int_tvol(MDIM,sm_NumBodies), int_rhog(MDIM,sm_NumBodies), int_xrhog(MDIM,sm_NumBodies), &
            int_rhou(MDIM,sm_NumBodies), int_xrhou(MDIM,sm_NumBodies) )

  annbody_type     = 0
  annbody_nparam   = 0
  annbody_param    = 0.
  r_bod = 0.; rd_bod = 0.; NwB_Nbod = 0.; TNB_bod = 0.;
  int_tvol = 0.; int_rhog = 0.; int_xrhog =0.; int_rhou = 0.; int_xrhou =0.;

  ! First populate body info
  do ibd = 1,sm_NumBodies
     if ( (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) .and. & 
          (sm_BodyInfo(ibd)%BodyType .eq. BODYTYPE_RIGID) ) then

        ! Parameters:
        annbody_type(ibd)     = sm_BodyInfo(ibd) % annbody_type
        annbody_nparam(ibd)   = sm_BodyInfo(ibd) % annbody_nparam
        if( annbody_nparam(ibd) .gt. 0 ) &
        annbody_param(1:annbody_nparam(ibd),ibd) = sm_BodyInfo(ibd) % annbody_param(1:annbody_nparam(ibd))

        ! Absolute positions, velocities and angular velocity:
        imaster = sm_BodyInfo(ibd) % borigin_node

        r_bod(IAXIS,ibd)=sm_BodyInfo(ibd)%x(imaster) + &
                         sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(IAXIS,imaster))
        r_bod(JAXIS,ibd)=sm_BodyInfo(ibd)%x(imaster) + &
                         sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(JAXIS,imaster))
#if NDIM == MDIM
        r_bod(KAXIS,ibd)=sm_BodyInfo(ibd)%x(imaster) + &
                         sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(KAXIS,imaster))
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


  call Timers_start("FluidVolIntRigid_allreduce")
  ! Now All reduce.... ughh expensive.
  ! Parameters
  aux_v = annbody_type
  call mpi_allreduce ( aux_v, annbody_type, sm_NumBodies, FLASH_INTEGER, &
                       MPI_SUM, sm_meshComm, ierr )
  aux_v = annbody_nparam
  call mpi_allreduce ( aux_v, annbody_nparam, sm_NumBodies, FLASH_INTEGER, &
                       MPI_SUM, sm_meshComm, ierr )

  ! Parameter list:
  aux_p = annbody_param
  call mpi_allreduce ( aux_p, annbody_param,RB_MAXANNPARAM*sm_NumBodies,FLASH_REAL, &
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
  call Timers_stop("FluidVolIntRigid_allreduce")


  ! Then, Block by block run over annalytical rigid bodies to 
  ! compute signed distance function on each velocity grid.
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Now do loop over Grid Blocks: 
  pt = 0
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get Blocks center coords:
     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get Blocks bsize:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     Vcell = product(del(1:NDIM))

#if NDIM == MDIM
     hcell = sqrt(del(IAXIS)**2. + del(JAXIS)**2. + del(KAXIS)**2.)
#else
     hcell = sqrt(del(IAXIS)**2. + del(JAXIS)**2. )
#endif

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
        zedge(k) = coord(KAXIS) - bsize(KAXIS)/2.0 + real(k - NGUARD -1)*del(KAXIS)
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
        yedge(j) = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j - NGUARD -1)*del(JAXIS)
        ycell(j) = yedge(j) + 0.5*del(JAXIS)
     enddo
     yedge(blkLimitsGC(HIGH,JAXIS)+1)=coord(JAXIS)-bsize(JAXIS)/2.0 + &
                                      real(blkLimitsGC(HIGH,JAXIS)-NGUARD)*del(JAXIS)

     ! X:
     do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
        xedge(i) = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i - NGUARD -1)*del(IAXIS)
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
           call Driver_abortFlash("Cannot Integrate sphere in 2D. Use RB_ANNDISC")
#endif     

        case ( RB_ANNELLIPSOID )

#if NDIM == MDIM              
           ! Ellipsoid:                    a                    b                    c       
           radius = max(annbody_param(1,ibd),annbody_param(2,ibd),annbody_param(3,ibd))/2.
           
           body_box(IAXIS,LOW)  = r_bod(IAXIS,ibd) - radius
           body_box(IAXIS,HIGH) = r_bod(IAXIS,ibd) + radius
           body_box(JAXIS,LOW)  = r_bod(JAXIS,ibd) - radius
           body_box(JAXIS,HIGH) = r_bod(JAXIS,ibd) + radius
           body_box(KAXIS,LOW)  = r_bod(KAXIS,ibd) - radius
           body_box(KAXIS,HIGH) = r_bod(KAXIS,ibd) + radius

#else
           call Driver_abortFlash("Cannot Integrate Ellipsoid in 2D.")
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

              ! Local position of lower cylinder center:  (p1) x----cen----x
              ! (p2)
              locp1(IAXIS,1) = -cyllen2

              ! Local position of higher cylinder center:
              locp2(IAXIS,1) =  cyllen2

           case(JAXIS) ! Cylinder direction in body axis yB  

              ! Local position of lower cylinder center:  (p1) x----cen----x
              ! (p2)
              locp1(JAXIS,1) = -cyllen2

              ! Local position of higher cylinder center:
              locp2(JAXIS,1) =  cyllen2

           case(KAXIS)

              ! Local position of lower cylinder center:  (p1) x----cen----x
              ! (p2)
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

        case ( RB_ANNRBC  )

           call Driver_abortFlash('RBC analytical geometry operations not coded yet.')

        case default

           call Driver_abortFlash('Analytical geometry type not known.')

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
        rhogx = 0.; rhogy = 0.; rhogz = 0.;
        rhou  = 0.; rhov  = 0.; rhow  = 0.;
        rhoDygx = 0.; rhoDzgx = 0.; 
        rhoDxgy = 0.; rhoDzgy = 0.;
        rhoDxgz = 0.; rhoDygz = 0.;
        rhoDyu = 0.; rhoDzu = 0.;
        rhoDxv = 0.; rhoDzv = 0.;
        rhoDxw = 0.; rhoDyw = 0.;
        Tvolx  = 0.; Tvoly  = 0.; Tvolz = 0.;        
        select case ( annbody_type(ibd) )
        case ( RB_ANNSPHERE )

#if NDIM == MDIM
           radius = annbody_param(1,ibd)/2.

           ! Z velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to Z veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zedge(1:MDIM,k)
                       vecnorm = sqrt(vec(IAXIS)**2. + vec(JAXIS)**2. + vec(KAXIS)**2.)

                       ! Inside sphere radius?
                       if ( vecnorm .lt. (radius+hcell) ) then

                          alphaijk=1.

                          ! Compute alphai 
                          if (vecnorm .ge. (radius - hcell)) then
                             ! Compute corners of Z face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xedgelow, yedgelow, zcelllow
                             xc(1:MDIM,1) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k-1)


                             ! xc2: xedgehigh, yedgelow, zcelllow
                             xc(1:MDIM,2) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k-1)


                             ! xc3: xedgelow, yedgehigh, zcelllow
                             xc(1:MDIM,3) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k-1)

                             ! xc4: xedgehigh, yedgehigh, zcelllow
                             xc(1:MDIM,4) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k-1)

                             ! xc5: xedgelow, yedgelow, zcellhigh
                             xc(1:MDIM,5) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k)


                             ! xc6: xedgehigh, yedgelow, zcellhigh
                             xc(1:MDIM,6) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k)


                             ! xc7: xedgelow, yedgehigh, zcellhigh
                             xc(1:MDIM,7) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k)
                          
                             ! xc8: xedgehigh, yedgehigh, zcellhigh
                             xc(1:MDIM,8) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k)

                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( xc(IAXIS,iph)**2.  + & 
                                                 xc(JAXIS,iph)**2.  + &
                                                 xc(KAXIS,iph)**2. )/radius -1.
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk =alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif


                          ! Total volume:
                          Tvolz = Tvolz + Vcell * alphaijk

                          ! rho*gz
                          rhogz= rhogz+ rhoi * sm_gravZ * alphaijk
                          
                          ! rho*w
                          rhow = rhow + rhoi * facezData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(y-ycm)*gz
                          rhoDygz= rhoDygz+ rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                    sm_gravZ * alphaijk

                          ! rho*(x-xcm)*gz
                          rhoDxgz= rhoDxgz+ rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                    sm_gravZ * alphaijk

                          ! rho*(y-ycm)*w
                          rhoDyw = rhoDyw + rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                    facezData(VELC_FACE_VAR,i,j,k) * alphaijk
                           
                          ! rho*(x-xcm)*w
                          rhoDxw = rhoDxw + rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                    facezData(VELC_FACE_VAR,i,j,k) * alphaijk  
                          if (k .ne. (blkLimits(HIGH,KAXIS)+1) ) pt(KAXIS) = pt(KAXIS) + 1
                       endif

                 enddo
              enddo
           enddo

           ! Y velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to Y veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_yedge(1:MDIM,j) + TBN_zcell(1:MDIM,k)
                       vecnorm = sqrt(vec(IAXIS)**2. + vec(JAXIS)**2. + vec(KAXIS)**2.)

                       ! Inside sphere radius?
                       if ( vecnorm .lt. (radius+hcell) ) then

                          alphaijk=1.

                          ! Compute alphai 
                          if (vecnorm .ge. (radius - hcell)) then
                             ! Compute corners of Y face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xedgelow, ycelllow, zedgelow
                             xc(1:MDIM,1) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc2: xedgehigh, ycelllow, zedgelow
                             xc(1:MDIM,2) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc3: xedgelow, ycellhigh, zedgelow
                             xc(1:MDIM,3) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc4: xedgehigh, ycellhigh, zedgelow
                             xc(1:MDIM,4) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc5: xedgelow, ycelllow, zedgehigh
                             xc(1:MDIM,5) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc6: xedgehigh, ycelllow, zedgehigh
                             xc(1:MDIM,6) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc7: xedgelow, ycellhigh, zedgehigh
                             xc(1:MDIM,7) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k+1)

                             ! xc8: xedgehigh, ycellhigh, zedgehigh
                             xc(1:MDIM,8) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( xc(IAXIS,iph)**2.  + &
                                                 xc(JAXIS,iph)**2.  + &
                                                 xc(KAXIS,iph)**2. )/radius -1.
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk =alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif

                          
                          ! Total volume:
                          Tvoly = Tvoly + Vcell * alphaijk

                          ! rho*gy
                          rhogy= rhogy+ rhoi * sm_gravY * alphaijk

                          ! rho*v
                          rhov = rhov + rhoi * faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(z-zcm)*gy
                          rhoDzgy= rhoDzgy+ rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   sm_gravY * alphaijk

                          ! rho*(x-xcm)*gy
                          rhoDxgy= rhoDxgy+ rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                   sm_gravY * alphaijk

                          ! rho*(z-zcm)*v
                          rhoDzv = rhoDzv + rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(x-xcm)*v
                          rhoDxv = rhoDxv + rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                   faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

                          if (j .ne. (blkLimits(HIGH,JAXIS)+1)) pt(JAXIS) = pt(JAXIS) + 1
                       endif

                 enddo
              enddo
           enddo

           ! X velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to X veloc point
                       vec(1:MDIM) = TBN_xedge(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zcell(1:MDIM,k)
                       vecnorm = sqrt(vec(IAXIS)**2. + vec(JAXIS)**2. + vec(KAXIS)**2.)

                       ! Inside sphere radius?
                       if ( vecnorm .lt. (radius+hcell) ) then

                          alphaijk=1.

                          ! Compute alphai 
                          if (vecnorm .ge. (radius - hcell)) then
                             ! Compute corners of Y face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xcelllow, yedgelow, zedgelow
                             xc(1:MDIM,1) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc2: xcellhigh, yedgelow, zedgelow
                             xc(1:MDIM,2) = TBN_xcell(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc3: xcelllow, yedgehigh, zedgelow
                             xc(1:MDIM,3) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc4: xcellhigh, yedgehigh, zedgelow
                             xc(1:MDIM,4) = TBN_xcell(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc5: xcelllow, yedgelow, zedgehigh
                             xc(1:MDIM,5) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc6: xcellhigh, yedgelow, zedgehigh
                             xc(1:MDIM,6) = TBN_xcell(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc7: xcelllow, yedgehigh, zedgehigh
                             xc(1:MDIM,7) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k+1)

                             ! xc8: xcellhigh, yedgehigh, zedgehigh
                             xc(1:MDIM,8) = TBN_xcell(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k+1)

                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( xc(IAXIS,iph)**2.  + &
                                                 xc(JAXIS,iph)**2.  + &
                                                 xc(KAXIS,iph)**2. )/radius -1.
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk=alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif

                          ! Total volume:
                          Tvolx = Tvolx + Vcell * alphaijk

                          ! rho*gx
                          rhogx= rhogx+ rhoi * sm_gravX * alphaijk

                          ! rho*u
                          rhou = rhou + rhoi * facexData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(z-zcm)*gx
                          rhoDzgx= rhoDzgx+ rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   sm_gravX * alphaijk

                          ! rho*(y-ycm)*gx
                          rhoDygx= rhoDygx+ rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                   sm_gravX * alphaijk

                          ! rho*(z-zcm)*u
                          rhoDzu = rhoDzu + rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   facexData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(y-ycm)*u
                          rhoDyu = rhoDyu + rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                   facexData(VELC_FACE_VAR,i,j,k) * alphaijk


                          if (i .ne. (blkLimits(HIGH,IAXIS)+1)) pt(IAXIS) = pt(IAXIS) + 1
                       endif

                 enddo
              enddo
           enddo
#else
           call Driver_abortFlash("Cannot Integrate sphere in 2D. Use RB_ANNDISC")

#endif 


        case ( RB_ANNELLIPSOID )

#if NDIM == MDIM

           radius1 = annbody_param(1,ibd)/2.
           radius2 = annbody_param(2,ibd)/2.
           radius3 = annbody_param(3,ibd)/2.
            

           ! Z velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to Z veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zedge(1:MDIM,k)

                       vecnorm = sqrt(   (vec(IAXIS)/(radius1+hcell))**2. &
                                       + (vec(JAXIS)/(radius2+hcell))**2. &
                                       + (vec(KAXIS)/(radius3+hcell))**2.)

                       ! Inside Ellipse?
                       if ( vecnorm .lt. 1. ) then

                          alphaijk=1.

                          vecnorm = sqrt(   (vec(IAXIS)/(radius1-hcell))**2. &
                                          + (vec(JAXIS)/(radius2-hcell))**2. &
                                          + (vec(KAXIS)/(radius3-hcell))**2.)                          
                          
                          ! Compute alphai 
                          if (vecnorm .ge. 1.) then
                             ! Compute corners of Z face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xedgelow, yedgelow, zcelllow
                             xc(1:MDIM,1) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k-1)


                             ! xc2: xedgehigh, yedgelow, zcelllow
                             xc(1:MDIM,2) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k-1)


                             ! xc3: xedgelow, yedgehigh, zcelllow
                             xc(1:MDIM,3) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k-1)

                             ! xc4: xedgehigh, yedgehigh, zcelllow
                             xc(1:MDIM,4) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k-1)

                             ! xc5: xedgelow, yedgelow, zcellhigh
                             xc(1:MDIM,5) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k)


                             ! xc6: xedgehigh, yedgelow, zcellhigh
                             xc(1:MDIM,6) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zcell(1:MDIM,k)


                             ! xc7: xedgelow, yedgehigh, zcellhigh
                             xc(1:MDIM,7) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k)
                          
                             ! xc8: xedgehigh, yedgehigh, zcellhigh
                             xc(1:MDIM,8) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zcell(1:MDIM,k)

                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( (xc(IAXIS,iph)/radius1)**2.  + & 
                                                 (xc(JAXIS,iph)/radius2)**2.  + &
                                                 (xc(KAXIS,iph)/radius3)**2. ) -1.
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk =alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif


                          ! Total volume:
                          Tvolz = Tvolz + Vcell * alphaijk

                          ! rho*gz
                          rhogz= rhogz+ rhoi * sm_gravZ * alphaijk
                          
                          ! rho*w
                          rhow = rhow + rhoi * facezData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(y-ycm)*gz
                          rhoDygz= rhoDygz+ rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                    sm_gravZ * alphaijk

                          ! rho*(x-xcm)*gz
                          rhoDxgz= rhoDxgz+ rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                    sm_gravZ * alphaijk

                          ! rho*(y-ycm)*w
                          rhoDyw = rhoDyw + rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                    facezData(VELC_FACE_VAR,i,j,k) * alphaijk
                           
                          ! rho*(x-xcm)*w
                          rhoDxw = rhoDxw + rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                    facezData(VELC_FACE_VAR,i,j,k) * alphaijk  
                          if (k .ne. (blkLimits(HIGH,KAXIS)+1) ) pt(KAXIS) = pt(KAXIS) + 1
                       endif

                 enddo
              enddo
           enddo

           ! Y velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to Y veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_yedge(1:MDIM,j) + TBN_zcell(1:MDIM,k)

                       vecnorm = sqrt(   (vec(IAXIS)/(radius1+hcell))**2. &
                                       + (vec(JAXIS)/(radius2+hcell))**2. &
                                       + (vec(KAXIS)/(radius3+hcell))**2.)

                       ! Inside Ellipse?
                       if ( vecnorm .lt. 1. ) then

                          alphaijk=1.

                          vecnorm = sqrt(   (vec(IAXIS)/(radius1-hcell))**2. &
                                          + (vec(JAXIS)/(radius2-hcell))**2. &
                                          + (vec(KAXIS)/(radius3-hcell))**2.)                          
                          
                          ! Compute alphai 
                          if (vecnorm .ge. 1.) then
                             ! Compute corners of Y face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xedgelow, ycelllow, zedgelow
                             xc(1:MDIM,1) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc2: xedgehigh, ycelllow, zedgelow
                             xc(1:MDIM,2) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc3: xedgelow, ycellhigh, zedgelow
                             xc(1:MDIM,3) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc4: xedgehigh, ycellhigh, zedgelow
                             xc(1:MDIM,4) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc5: xedgelow, ycelllow, zedgehigh
                             xc(1:MDIM,5) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc6: xedgehigh, ycelllow, zedgehigh
                             xc(1:MDIM,6) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc7: xedgelow, ycellhigh, zedgehigh
                             xc(1:MDIM,7) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k+1)

                             ! xc8: xedgehigh, ycellhigh, zedgehigh
                             xc(1:MDIM,8) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( (xc(IAXIS,iph)/radius1)**2.  + & 
                                                 (xc(JAXIS,iph)/radius2)**2.  + &
                                                 (xc(KAXIS,iph)/radius3)**2. ) -1.
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk =alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif

                          
                          ! Total volume:
                          Tvoly = Tvoly + Vcell * alphaijk

                          ! rho*gy
                          rhogy= rhogy+ rhoi * sm_gravY * alphaijk

                          ! rho*v
                          rhov = rhov + rhoi * faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(z-zcm)*gy
                          rhoDzgy= rhoDzgy+ rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   sm_gravY * alphaijk

                          ! rho*(x-xcm)*gy
                          rhoDxgy= rhoDxgy+ rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                   sm_gravY * alphaijk

                          ! rho*(z-zcm)*v
                          rhoDzv = rhoDzv + rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(x-xcm)*v
                          rhoDxv = rhoDxv + rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                                   faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

                          if (j .ne. (blkLimits(HIGH,JAXIS)+1)) pt(JAXIS) = pt(JAXIS) + 1
                       endif

                 enddo
              enddo
           enddo

           ! X velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to X veloc point
                       vec(1:MDIM) = TBN_xedge(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zcell(1:MDIM,k)

                       vecnorm = sqrt(   (vec(IAXIS)/(radius1+hcell))**2. &
                                       + (vec(JAXIS)/(radius2+hcell))**2. &
                                       + (vec(KAXIS)/(radius3+hcell))**2.)

                       ! Inside Ellipse?
                       if ( vecnorm .lt. 1. ) then

                          alphaijk=1.

                          vecnorm = sqrt(   (vec(IAXIS)/(radius1-hcell))**2. &
                                          + (vec(JAXIS)/(radius2-hcell))**2. &
                                          + (vec(KAXIS)/(radius3-hcell))**2.)                          
                          
                          ! Compute alphai 
                          if (vecnorm .ge. 1.) then
                             ! Compute corners of Y face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xcelllow, yedgelow, zedgelow
                             xc(1:MDIM,1) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc2: xcellhigh, yedgelow, zedgelow
                             xc(1:MDIM,2) = TBN_xcell(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc3: xcelllow, yedgehigh, zedgelow
                             xc(1:MDIM,3) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc4: xcellhigh, yedgehigh, zedgelow
                             xc(1:MDIM,4) = TBN_xcell(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc5: xcelllow, yedgelow, zedgehigh
                             xc(1:MDIM,5) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc6: xcellhigh, yedgelow, zedgehigh
                             xc(1:MDIM,6) = TBN_xcell(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k+1)


                             ! xc7: xcelllow, yedgehigh, zedgehigh
                             xc(1:MDIM,7) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k+1)

                             ! xc8: xcellhigh, yedgehigh, zedgehigh
                             xc(1:MDIM,8) = TBN_xcell(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k+1)

                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( (xc(IAXIS,iph)/radius1)**2.  + & 
                                                 (xc(JAXIS,iph)/radius2)**2.  + &
                                                 (xc(KAXIS,iph)/radius3)**2. ) -1.
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk=alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif

                          ! Total volume:
                          Tvolx = Tvolx + Vcell * alphaijk

                          ! rho*gx
                          rhogx= rhogx+ rhoi * sm_gravX * alphaijk

                          ! rho*u
                          rhou = rhou + rhoi * facexData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(z-zcm)*gx
                          rhoDzgx= rhoDzgx+ rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   sm_gravX * alphaijk

                          ! rho*(y-ycm)*gx
                          rhoDygx= rhoDygx+ rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                   sm_gravX * alphaijk

                          ! rho*(z-zcm)*u
                          rhoDzu = rhoDzu + rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                                   facexData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(y-ycm)*u
                          rhoDyu = rhoDyu + rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                                   facexData(VELC_FACE_VAR,i,j,k) * alphaijk


                          if (i .ne. (blkLimits(HIGH,IAXIS)+1)) pt(IAXIS) = pt(IAXIS) + 1
                       endif

                 enddo
              enddo
           enddo
#else
           call Driver_abortFlash("Cannot Integrate Ellipsoid in 2D.")

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

#if NDIM == MDIM
           ! Z velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to Z veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zedge(1:MDIM,k)

                       inlenflg= abs(vec(axis_1)) .lt. (cyllen2 + hcell)
                       vecnorm = sqrt(vec(axis_2)**2. + vec(axis_3)**2.)

                       ! Inside cylinder?
                       if ( inlenflg .and. (vecnorm .lt. (radius + hcell)) ) then

                          ! Total volume:
                          Tvolz = Tvolz + Vcell

                          ! rho*gz
                          rhogz= rhogz+ rhoi * sm_gravZ

                          ! rho*w
                          rhow = rhow + rhoi * facezData(VELC_FACE_VAR,i,j,k)

                          ! rho*(y-ycm)*gz
                          rhoDygz= rhoDygz+ rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * sm_gravZ

                          ! rho*(x-xcm)*gz
                          rhoDxgz= rhoDxgz+ rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * sm_gravZ

                          ! rho*(y-ycm)*w
                          rhoDyw = rhoDyw + rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * facezData(VELC_FACE_VAR,i,j,k)

                          ! rho*(x-xcm)*w
                          rhoDxw = rhoDxw + rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * facezData(VELC_FACE_VAR,i,j,k) 

                          if (k .lt. (blkLimits(HIGH,KAXIS)+1)) pt(KAXIS) = pt(KAXIS) + 1
                       endif

                 enddo
              enddo
           enddo
#endif
           ! Y velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to Y veloc point
                       vec(1:MDIM) = TBN_xcell(1:MDIM,i) + TBN_yedge(1:MDIM,j) + TBN_zcell(1:MDIM,k)

                       inlenflg= abs(vec(axis_1)) .lt. (cyllen2 + hcell)
                       vecnorm = sqrt(vec(axis_2)**2. + vec(axis_3)**2.)

                       ! Inside cylinder?
                       if ( inlenflg .and. (vecnorm .lt. (radius + hcell)) ) then

                          alphaijk = 1.
#if NDIM == MDIM

#else
                          ! Compute alphai 
                          if (vecnorm .ge. (radius - hcell)) then
                             ! Compute corners of Y face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xedgelow, ycelllow, zedgelow
                             xc(1:MDIM,1) = TBN_xedge(1:MDIM,i)   + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k)

                                    
                             ! xc2: xedgehigh, ycelllow, zedgelow
                             xc(1:MDIM,2) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j-1) + &
                                            TBN_zedge(1:MDIM,k)
  

                             ! xc3: xedgelow, ycellhigh, zedgelow
                             xc(1:MDIM,3) = TBN_xedge(1:MDIM,i  ) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc4: xedgehigh, ycellhigh, zedgelow
                             xc(1:MDIM,4) = TBN_xedge(1:MDIM,i+1) + &
                                            TBN_ycell(1:MDIM,j  ) + &
                                            TBN_zedge(1:MDIM,k)


                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( xc(axis_2,iph)**2.  + &
                                                 xc(axis_3,iph)**2. )/radius - 1. 
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk = alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif

#endif


                          ! Total volume:
                          Tvoly = Tvoly + Vcell * alphaijk

                          ! rho*gy
                          rhogy= rhogy+ rhoi * sm_gravY * alphaijk

                          ! rho*v
                          rhov = rhov + rhoi * faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(x-xcm)*gy
                          rhoDxgy= rhoDxgy+ rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * sm_gravY * alphaijk

                          ! rho*(x-xcm)*v
                          rhoDxv = rhoDxv + rhoi * (xcell(i)-r_bod(IAXIS,ibd)) * &
                                            faceyData(VELC_FACE_VAR,i,j,k) * alphaijk

#if NDIM == MDIM
                          ! rho*(z-zcm)*gy            
                          rhoDzgy= rhoDzgy+ rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * sm_gravY * alphaijk

                          ! rho*(z-zcm)*v
                          rhoDzv = rhoDzv + rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                            faceyData(VELC_FACE_VAR,i,j,k) * alphaijk
#endif
                          if (j .lt. (blkLimits(HIGH,JAXIS)+1)) pt(JAXIS) = pt(JAXIS) + 1

                       endif

                 enddo
              enddo
           enddo
           ! X velocity:
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                       ! Distance from body center of mass to X veloc point
                       vec(1:MDIM) = TBN_xedge(1:MDIM,i) + TBN_ycell(1:MDIM,j) + TBN_zcell(1:MDIM,k)

                       inlenflg= abs(vec(axis_1)) .lt. (cyllen2 + hcell)
                       vecnorm = sqrt(vec(axis_2)**2. + vec(axis_3)**2.)

                       ! Inside cylinder?
                       if ( inlenflg .and. (vecnorm .lt. (radius + hcell)) ) then

                          alphaijk = 1.
#if NDIM == MDIM

#else
                          ! Compute alphai 
                          if (vecnorm .ge. (radius - hcell)) then
                             ! Compute corners of Y face control volume:
                             ! if edge: low = i, high = i+1
                             !    cell: low = i-1, high = i
                             ! xc1: xcelllow, yedgelow, zedgelow
                             xc(1:MDIM,1) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc2: xcellhigh, yedgelow, zedgelow
                             xc(1:MDIM,2) = TBN_xcell(1:MDIM,i)   + &
                                            TBN_yedge(1:MDIM,j)   + &
                                            TBN_zedge(1:MDIM,k)


                             ! xc3: xcelllow, yedgehigh, zedgelow
                             xc(1:MDIM,3) = TBN_xcell(1:MDIM,i-1) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k)

                             ! xc4: xcellhigh, yedgehigh, zedgelow
                             xc(1:MDIM,4) = TBN_xcell(1:MDIM,i  ) + &
                                            TBN_yedge(1:MDIM,j+1) + &
                                            TBN_zedge(1:MDIM,k)


                             ! Signed distance function: -ve inside body:
                             do iph = 1,2**NDIM
                                phi(iph) = sqrt( xc(axis_2,iph)**2.  + &
                                                 xc(axis_3,iph)**2. )/radius - 1.
                             enddo

                             ! Solid Volume fraction: Kempe & Frohlich JCP 2012
                             alphaijk = 0.
                             do iph = 1,2**NDIM
                                alphaijk = alphaijk-phi(iph)*0.5*(sign(1.,-phi(iph))+1.)
                             enddo
                             alphaijk = alphaijk/sum(abs(phi(1:2**NDIM)))

                          endif

#endif


                          ! Total volume:
                          Tvolx = Tvolx + Vcell * alphaijk

                          ! rho*gx
                          rhogx= rhogx+ rhoi * sm_gravX * alphaijk

                          ! rho*u
                          rhou = rhou + rhoi * facexData(VELC_FACE_VAR,i,j,k) * alphaijk

                          ! rho*(y-ycm)*gx
                          rhoDygx= rhoDygx+ rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * sm_gravX * alphaijk

                          ! rho*(y-ycm)*u
                          rhoDyu = rhoDyu + rhoi * (ycell(j)-r_bod(JAXIS,ibd)) * &
                                            facexData(VELC_FACE_VAR,i,j,k) * alphaijk

#if NDIM == MDIM
                          ! rho*(z-zcm)*gx
                          rhoDzgx= rhoDzgx+ rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * sm_gravX * alphaijk

                          ! rho*(z-zcm)*u
                          rhoDzu = rhoDzu + rhoi * (zcell(k)-r_bod(KAXIS,ibd)) * &
                                            facexData(VELC_FACE_VAR,i,j,k) * alphaijk
#endif
                          if (i .lt. (blkLimits(HIGH,IAXIS)+1)) pt(IAXIS) = pt(IAXIS) + 1

                       endif

                 enddo
              enddo
           enddo

        case ( RB_ANNRBC  )

           call Driver_abortFlash('RBC analytical geometry operations not coded yet.')

        case default

           call Driver_abortFlash('Analytical geometry type not known.')

        end select

        ! Here add to body Integrals:
        ! ...

        int_tvol(IAXIS,ibd)  =  int_tvol(IAXIS,ibd) + Tvolx
        int_tvol(JAXIS,ibd)  =  int_tvol(JAXIS,ibd) + Tvoly

        int_rhog(IAXIS,ibd)  =  int_rhog(IAXIS,ibd) + rhogx*Vcell
        int_rhog(JAXIS,ibd)  =  int_rhog(JAXIS,ibd) + rhogy*Vcell

        int_xrhog(KAXIS,ibd) = int_xrhog(KAXIS,ibd) + (rhoDxgy-rhoDygx)*Vcell

        int_rhou(IAXIS,ibd)  =  int_rhou(IAXIS,ibd) + rhou*Vcell
        int_rhou(JAXIS,ibd)  =  int_rhou(JAXIS,ibd) + rhov*Vcell

        int_xrhou(KAXIS,ibd) = int_xrhou(KAXIS,ibd) + (rhoDxv-rhoDyu)*Vcell

#if NDIM == MDIM
        int_tvol(KAXIS,ibd)  =  int_tvol(KAXIS,ibd) + Tvolz

        int_rhog(KAXIS,ibd)  =  int_rhog(KAXIS,ibd) + rhogz*Vcell
        int_xrhog(IAXIS,ibd) = int_xrhog(IAXIS,ibd) + (rhoDygz-rhoDzgy)*Vcell
        int_xrhog(JAXIS,ibd) = int_xrhog(JAXIS,ibd) + (rhoDzgx-rhoDxgz)*Vcell

        int_rhou(KAXIS,ibd)  =  int_rhou(KAXIS,ibd) + rhow*Vcell
        int_xrhou(IAXIS,ibd) = int_xrhou(IAXIS,ibd) + (rhoDyw-rhoDzv)*Vcell
        int_xrhou(JAXIS,ibd) = int_xrhou(JAXIS,ibd) + (rhoDzu-rhoDxw)*Vcell
#endif


     enddo ! Loop Bodies

     ! Release Pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo ! Loop Blocks

  !write(*,*) 'rhoU=',sm_MeshMe,int_rhou(IAXIS:JAXIS,1)

  ! All Reduce integrals:
  aux_r = int_tvol
  call mpi_allreduce ( aux_r, int_tvol, MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr )

  aux_r = int_rhog
  call mpi_allreduce ( aux_r, int_rhog, MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr )

  aux_r = int_xrhog
  call mpi_allreduce ( aux_r, int_xrhog, MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr )
  
  aux_r = int_rhou
  call mpi_allreduce ( aux_r, int_rhou, MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr )

  aux_r = int_xrhou
  call mpi_allreduce ( aux_r, int_xrhou, MDIM*sm_NumBodies, FLASH_REAL, &
                       MPI_SUM, sm_meshComm, ierr )


  ! Now sums of integrals go back to Body Masters:
  do ibd = 1,sm_NumBodies
     if ( (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) .and. &
          (sm_BodyInfo(ibd)%BodyType .eq. BODYTYPE_RIGID) ) then

        ! Compute Analytical volume:
        select case ( annbody_type(ibd) )

        case ( RB_ANNSPHERE )

#if NDIM == MDIM              
           ! Sphere:
           radius   = annbody_param(1,ibd)/2.
           body_vol = 4./3.* PI * radius**3.
#else
           call Driver_abortFlash("Cannot Integrate sphere in 2D. Use RB_ANNDISC")
#endif           

        case ( RB_ANNELLIPSOID )

#if NDIM == MDIM              
           ! Ellipsoid:
           radius1   = annbody_param(1,ibd)/2.
           radius2   = annbody_param(2,ibd)/2.
           radius3   = annbody_param(3,ibd)/2.
           
           body_vol = 4./3.* PI * radius1 * radius2 * radius3
#else
           call Driver_abortFlash("Cannot Integrate Ellipsoid in 2D.")
#endif     

        case ( RB_ANNDISC )

           ! Cylinder
           radius = annbody_param(2,ibd)/2.
#if NDIM == MDIM
           cyllen = annbody_param(3,ibd)
#else
           cyllen = 1.
#endif            
           body_vol = PI * radius**2. * cyllen
 
        case ( RB_ANNRBC  )
           call Driver_abortFlash('RBC analytical geometry operations not coded yet.')
        case default
           call Driver_abortFlash('Analytical geometry type not known.')
        end select

        ! Volume factor:
        vol_fctx = body_vol/int_tvol(IAXIS,ibd)
        vol_fcty = body_vol/int_tvol(JAXIS,ibd)
      
        !write(*,*) 'Volumes X,Y,ANN=',int_tvol(IAXIS,ibd),int_tvol(JAXIS,ibd),int_tvol(JAXIS,ibd)

        ! Assign to BodyInfo data:
        sm_BodyInfo(ibd)%int_rhog(IAXIS) = int_rhog(IAXIS,ibd)*vol_fctx
        sm_BodyInfo(ibd)%int_rhog(JAXIS) = int_rhog(JAXIS,ibd)*vol_fcty


        sm_BodyInfo(ibd)%int_rhou(IAXIS,CONSTANT_ONE) = int_rhou(IAXIS,ibd)!*vol_fctx
        sm_BodyInfo(ibd)%int_rhou(JAXIS,CONSTANT_ONE) = int_rhou(JAXIS,ibd)!*vol_fcty

#if NDIM == MDIM 
        vol_fctz = body_vol/int_tvol(KAXIS,ibd)

        sm_BodyInfo(ibd)%int_rhog(KAXIS) = int_rhog(KAXIS,ibd)*vol_fctz
        sm_BodyInfo(ibd)%int_rhou(KAXIS,CONSTANT_ONE) = int_rhou(KAXIS,ibd)!*vol_fctz

        !write(*,*) 'Volumes X,Y,Z,ANN=',int_tvol(IAXIS,ibd),int_tvol(JAXIS,ibd), &
        !                                int_tvol(KAXIS,ibd),body_vol
        !write(*,*) 'Volume rhou rhog Z=',int_rhou(KAXIS,ibd),int_rhog(KAXIS,ibd)*vol_fctz

#endif         

        sm_BodyInfo(ibd)%int_xrhog(1:MDIM)= int_xrhog(1:MDIM,ibd)*(vol_fctx+vol_fcty+vol_fctz)/3.
        sm_BodyInfo(ibd)%int_xrhou(1:MDIM,CONSTANT_ONE)= int_xrhou(1:MDIM,ibd)

     endif
  enddo  


  ! Total Points in u,v,w integrated:
  auxpt = pt
  call mpi_allreduce( auxpt, pt, MDIM, FLASH_INTEGER, MPI_SUM, sm_meshComm, ierr)
  if (sm_meshMe .eq. MASTER_PE) write(*,*) 'Total Internal Points Integrated u,v,w=',pt(IAXIS:KAXIS)


  ! Deallocate:
  deallocate( annbody_type , annbody_nparam )
  deallocate( annbody_param )
  deallocate( r_bod , rd_bod , NwB_Nbod )
  deallocate( TNB_bod )
  deallocate( aux_v , aux_p , aux_r , aux_t )

  deallocate( int_rhog, int_xrhog, int_rhou, int_xrhou )


  return
end subroutine sm_assemble_FluidVolInt_rigid
