

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

#define TWO_POINTSP 1
#define NORMAL_GRAD_CORR
!#define TEST_COMPARE

subroutine ib_distributedForces(blockID, particleData, vortx, vorty, vortz)

  use Grid_Data, only : gr_meshMe

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,      &
                             Grid_getDeltas, Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize

  use ImBound_Data, only : ib_nu,ib_stencil,ib_alphax,ib_alphay,ib_alphaz,ib_dt

  use ib_interface, only : ib_stencils,ib_getInterpFunc

  use Driver_interface, only : Driver_abortFlash

  use Grid_data,ONLY : gr_imin,gr_jmin,gr_kmin

  use IncompNS_data, ONLY : ins_gravX,ins_gravY,ins_gravZ

  implicit none
  integer, intent(IN) :: blockID
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortx 
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vorty
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortz
  real, intent(INOUT) :: particleData(NPART_PROPS)

  ! Local Variables....
  real :: xp,yp,zp,zL,h,hl,dx,dy,dz,dsx,dsy,dsz,ubd,vbd,wbd,ubdd,vbdd,wbdd,nxp,nyp,nzp
  real, dimension(MDIM) :: xbe,del,coord,bsize,np

  integer, parameter, dimension(MDIM)      :: grdip = (/ CENTER, CENTER, CENTER /)
  real, parameter, dimension(MDIM)      :: dlip = (/ 0.5, 0.5, 0.5 /)


#ifdef TANGENT_WITH_VORTICITY
  integer, parameter, dimension(MDIM,MDIM) :: grdnu = &
           RESHAPE( (/  FACES, FACES,CENTER, CENTER, FACES, FACES, FACES,CENTER, FACES /), (/MDIM,MDIM /))
                    !    k (wz)               i (wx)                j (wy)
  real, parameter, dimension(MDIM,MDIM) :: dlnu = &
           RESHAPE( (/ 0., 0., 0.5, 0.5, 0., 0., 0., 0.5, 0. /), (/MDIM,MDIM /))  
                    !   k             i           j
#else
  integer, parameter, dimension(MDIM,MDIM) :: grdu = &
           RESHAPE( (/  FACES,CENTER,CENTER, CENTER, FACES,CENTER, CENTER,CENTER, FACES /), (/MDIM,MDIM /))
                    !    u                   v                     w
  real, parameter, dimension(MDIM,MDIM) :: dlu = &
           RESHAPE( (/ 0., 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0. /), (/MDIM,MDIM /))
                    !  u             v             w
  real :: ui,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,exx,eyy,ezz,exy,exz,eyz
  real :: ue,ve,we,ues,ves,wes,ven,ups,vps,wps,vpn,dun,dvn,dwn,dune,dvne,dwne

  real :: normt,tx,ty,tz,dpdxt,dpdx,dpdy,dpdz

  real :: zpres2,xbe2(MDIM),zpres3,xbe3(MDIM)
#endif

  integer :: presflag,gridfl(MDIM),i,idim,gridind,nkij
  integer :: ielem(ib_stencil,MDIM,CONSTANT_TWO)
  real :: dpdn,eps
  real :: delaux(MDIM),xyz_stencil(ib_stencil,MDIM),phile(ib_stencil,NDIM+1)

  real :: zpres,p_i,zv(MDIM),nuwx,nuwy,nuwz

  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData

#ifdef TANGENT_WITH_VORTICITY
  integer, parameter :: derivflag = 0 ! Give interpolation functions.
#else
  integer, parameter :: derivflag = 1 ! Give interpolation functions and their derivatives
#endif

#ifdef NORMAL_GRAD_CORR
  real :: alpha
  real :: imp(3,3,3)
  real :: ddvel(3,3,3), ddp(3,3)
  real :: tp(3)
  real :: ddpdxdx, ddpdxdy, ddpdxdz, ddpdydy, ddpdydz, ddpdzdz
  real :: ddudxdx, ddudxdy, ddudxdz, ddudydy, ddudydz, ddudzdz
  real :: ddvdxdx, ddvdxdy, ddvdxdz, ddvdydy, ddvdydz, ddvdzdz
  real :: ddwdxdx, ddwdxdy, ddwdxdz, ddwdydy, ddwdydz, ddwdzdz
  real :: dpdxn, zL2
  real :: fvx_l, fvy_l, fvz_l
  real :: fvx_c, fvy_c, fvz_c, fv_c
  real :: fvx_cc, fvy_cc, fvz_cc, fv_cc
  real :: w_l(3), w_c(3), w_cc(3)
  integer :: ii, jj, kk
#endif

  real ::nu, dt

  nu = ib_nu
  dt = ib_dt

#ifndef TANGENT_WITH_VORTICITY  ! Shizhao Wang
  ue =0.; ve=0.; we=0.; 

  dudx=0.; dudy=0.; dudz=0.;
  dvdx=0.; dvdy=0.; dvdz=0.; 
  dwdx=0.; dwdy=0.; dwdz=0.;

  exx=0;  exy=0.; exz=0.;
  eyy=0.; eyz=0.;
  ezz=0.;
#endif

  ! Get dx,dy
  call Grid_getDeltas(blockID,del)
  call Grid_getBlkCenterCoords(blockID,coord)
  call Grid_getBlkPhysicalSize(blockID,bsize)

  dx=Del(IAXIS)
  dy=Del(JAXIS)
#if NDIM == 3
  dz=Del(KAXIS)
#else
  dz=1.
#endif

  ! Particle data:
  xp  = particleData(POSX_PART_PROP)
  yp  = particleData(POSY_PART_PROP)
  ubd = particleData(VELX_PART_PROP)
  vbd = particleData(VELY_PART_PROP)
  ubdd= particleData(ACCX_PART_PROP)
  vbdd= particleData(ACCY_PART_PROP)
  nxp = particleData(NMLX_PART_PROP)
  nyp = particleData(NMLY_PART_PROP)

  np(IAXIS) = nxp
  np(JAXIS) = nyp

#if NDIM == 3
  zp  = particleData(POSZ_PART_PROP)
  wbd = particleData(VELZ_PART_PROP)
  wbdd= particleData(ACCZ_PART_PROP)
  nzp = particleData(NMLZ_PART_PROP)
#else
  zp  = 0.
  wbd = 0.
  wbdd= 0.
  nzp = 0.
#endif
  np(KAXIS) = nzp


  ! Function Gimme h:
  ! call ib_normaldistance(dx,dy,dz,h)
  dsx = ib_alphax*dx
  dsy = ib_alphay*dy
#if NDIM == 2
  !h   = 0.5*(dsx+dsy) 
  h   = 1.*sqrt( (dsx*nxp)**2. + (dsy*nyp)**2. )
  eps = 1.E-2*MIN(dx,dy)  
#elif NDIM == 3
  dsz = ib_alphaz*dz
  !h   = 1./3.*(dsx+dsy+dsz)
  h   = 1.*sqrt( (dsx*nxp)**2. + (dsy*nyp)**2. + (dsz*nzp)**2. )
  eps = 1.E-2*MIN(dx,dy,dz)
#endif

             
  ! External Point Position:
  xbe(IAXIS) = xp + nxp*h
  xbe(JAXIS) = yp + nyp*h
#if NDIM == 3
  xbe(KAXIS) = zp + nzp*h
#else
  xbe(KAXIS) = 0.
#endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Add by Shizhao Wang, Oct 22, 2014
! To take account the variation of pressure in normal idrection on coarse mesh
#ifdef NORMAL_GRAD_CORR
  
  ! The normal external probe point is set to be 2.0dx from the wall
  alpha = real(NDIM)

  ddvel = 0.
  ddp = 0.
  tp = 0.

   ! Initial the flow in z direction. This is necessary for the 2D case
   we   = 0.; dwdx = 0.; dwdy =0.;
   dwdz = 0.
   ddwdxdx = 0.
   ddwdxdy = 0.
   ddwdxdz = 0.
   ddwdydy = 0.
   ddwdydz = 0.
   ddwdzdz = 0.
   ddpdxdz = 0.; ddpdydz = 0.; ddpdzdz =0.;
   ddudxdz = 0.; ddudydz = 0.; ddudzdz =0.;
   ddvdxdz = 0.; ddvdydz = 0.; ddvdzdz =0.;

  ! Function Gimme h:
  ! call ib_normaldistance(dx,dy,dz,h)
  dsx = alpha*dx
  dsy = alpha*dy
#if NDIM == 2
  !h   = 0.5*(dsx+dsy) 
  h   = 1.*sqrt( (dsx*nxp)**2. + (dsy*nyp)**2. )
  eps = 1.E-2*MIN(dx,dy)  
#elif NDIM == 3
  dsz = alpha*dz
  !h   = 1./3.*(dsx+dsy+dsz)
  h   = 1.*sqrt( (dsx*nxp)**2. + (dsy*nyp)**2. + (dsz*nzp)**2. )
  eps = 1.E-2*MIN(dx,dy,dz)
#endif
             
  ! External Point Position:
  xbe(IAXIS) = xp + nxp*h
  xbe(JAXIS) = yp + nyp*h
#if NDIM == 3
  xbe(KAXIS) = zp + nzp*h
#else
  xbe(KAXIS) = 0.
#endif

  zpres = 0.
  zv(1:MDIM) = 0.
  zL= 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
#ifdef TANGENT_WITH_VORTICITY
!     nkij = 1 + (1-presflag)*(2*NDIM-MDIM-1)
#else
     nkij = 1 + (1-presflag)*(NDIM-1) ! in 2D d/dx,d/dy; in 3D d/dx,d/dy,d/dz
#endif

     do gridind = 1,nkij

#ifdef TANGENT_WITH_VORTICITY
        ! Define Grids in case of vorticity and pressure:
!        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdnu(:,gridind)
        ! Auxiliary deltas
!        do idim = 1,NDIM
!           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlnu(idim,gridind))*del(idim)
!        enddo
#else
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlu(idim,gridind))*del(idim)
        enddo
#endif

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe,np,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe,xyz_stencil,del,derivflag,phile)
    
        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

#ifndef TANGENT_WITH_VORTICITY
           dpdx = 0.; dpdy =0.;
           dpdz = 0.
           ddpdxdx = 0.
           ddpdxdy = 0.
           ddpdxdz = 0.
           ddpdydy = 0.
           ddpdydz = 0.
           ddpdzdz = 0. 
#endif

           do i = 1 , ib_stencil      
              imp(1,1,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                          ielem(i,JAXIS,presflag+1)-1, &
                                          ielem(i,KAXIS,presflag+1));
              imp(2,1,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                        ielem(i,JAXIS,presflag+1)-1, &
                                        ielem(i,KAXIS,presflag+1));
              imp(3,1,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                          ielem(i,JAXIS,presflag+1)-1, &
                                          ielem(i,KAXIS,presflag+1));
 
              imp(1,2,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                        ielem(i,JAXIS,presflag+1), &
                                        ielem(i,KAXIS,presflag+1));
              imp(2,2,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1));  
              imp(3,2,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                        ielem(i,JAXIS,presflag+1), &
                                        ielem(i,KAXIS,presflag+1));

              imp(1,3,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                          ielem(i,JAXIS,presflag+1)+1, &
                                          ielem(i,KAXIS,presflag+1));
              imp(2,3,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                        ielem(i,JAXIS,presflag+1)+1, &
                                        ielem(i,KAXIS,presflag+1));
              imp(3,3,2) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                          ielem(i,JAXIS,presflag+1)+1, &
                                          ielem(i,KAXIS,presflag+1));

               zpres = zpres + phile(i,1)*imp(2,2,2)
               dpdx  = dpdx  + phile(i,1)*(imp(3,2,2)-imp(1,2,2))/(2.*dx)
               dpdy  = dpdy  + phile(i,1)*(imp(2,3,2)-imp(2,1,2))/(2.*dy)
               ddpdxdx = ddpdxdx + phile(i,1)*(imp(3,2,2)-2.*imp(2,2,2)+imp(1,2,2))/(dx*dx)
               ddpdydy = ddpdydy + phile(i,1)*(imp(2,3,2)-2.*imp(2,2,2)+imp(2,1,2))/(dy*dy)
               ddpdxdy = ddpdxdy + phile(i,1)*(imp(3,3,2)-imp(3,1,2)-imp(1,3,2)+imp(1,1,2))/(4.*dx*dy)

#if NDIM == MDIM
!              imp(1,1,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                          ielem(i,JAXIS,presflag+1)-1, &
!                                          ielem(i,KAXIS,presflag+1)-1);
              imp(2,1,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                        ielem(i,JAXIS,presflag+1)-1, &
                                        ielem(i,KAXIS,presflag+1)-1);
!              imp(3,1,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                          ielem(i,JAXIS,presflag+1)-1, &
!                                          ielem(i,KAXIS,presflag+1)-1);
 
              imp(1,2,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                        ielem(i,JAXIS,presflag+1), &
                                        ielem(i,KAXIS,presflag+1)-1);
              imp(2,2,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1)-1);  
              imp(3,2,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                        ielem(i,JAXIS,presflag+1), &
                                        ielem(i,KAXIS,presflag+1)-1);

!              imp(1,3,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                          ielem(i,JAXIS,presflag+1)+1, &
!                                          ielem(i,KAXIS,presflag+1)-1);
              imp(2,3,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                        ielem(i,JAXIS,presflag+1)+1, &
                                        ielem(i,KAXIS,presflag+1)-1);
!              imp(3,3,1) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                          ielem(i,JAXIS,presflag+1)+1, &
!                                          ielem(i,KAXIS,presflag+1)-1);
              !------
!              imp(1,1,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                          ielem(i,JAXIS,presflag+1)-1, &
!                                          ielem(i,KAXIS,presflag+1)+1);
              imp(2,1,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                        ielem(i,JAXIS,presflag+1)-1, &
                                        ielem(i,KAXIS,presflag+1)+1);
!              imp(3,1,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                          ielem(i,JAXIS,presflag+1)-1, &
!                                          ielem(i,KAXIS,presflag+1)+1);
 
              imp(1,2,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                        ielem(i,JAXIS,presflag+1), &
                                        ielem(i,KAXIS,presflag+1)+1);
              imp(2,2,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1)+1);  
              imp(3,2,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                        ielem(i,JAXIS,presflag+1), &
                                        ielem(i,KAXIS,presflag+1)+1);

!              imp(1,3,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                          ielem(i,JAXIS,presflag+1)+1, &
!                                          ielem(i,KAXIS,presflag+1)+1);
              imp(2,3,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                        ielem(i,JAXIS,presflag+1)+1, &
                                        ielem(i,KAXIS,presflag+1)+1);
!              imp(3,3,3) = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                          ielem(i,JAXIS,presflag+1)+1, &
!                                          ielem(i,KAXIS,presflag+1)+1);

               dpdz  = dpdz  + phile(i,1)*(imp(2,2,3)-imp(2,2,1))/(2.*dz)
               ddpdzdz = ddpdzdz + phile(i,1)*(imp(2,2,3)-2.*imp(2,2,2)+imp(2,2,1))/(dz*dz)
               ddpdxdz = ddpdxdz + phile(i,1)*(imp(3,2,3)-imp(3,2,1)-imp(1,2,3)+imp(1,2,1))/(4.*dx*dz)
               ddpdydz = ddpdydz + phile(i,1)*(imp(2,3,3)-imp(2,3,1)-imp(2,1,3)+imp(2,1,1))/(4.*dy*dz)

#endif             
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           ! Get Pressure approximation at surface marker: Acceleration +
           ! gravity effects.
           dpdn = -(     ubdd*nxp +      vbdd*nyp +      wbdd*nzp) + & ! -rho*Du/Dt * n 
                   (ins_gravX*nxp + ins_gravY*nyp + ins_gravZ*nzp);    ! +rho*    g * n
           zL = zpres - dpdn*h;

        else                                         ! Tangent stress

#ifdef TANGENT_WITH_VORTICITY
!           select case (gridind)
#else

           select case(gridind)
           case(1)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facexData,FACEX)

           ! U velocity derivatives on external point:
           ue = 0.; dudx = 0.; dudy =0.;
           dudz = 0.
           ddudxdx = 0.
           ddudxdy = 0.
           ddudxdz = 0.
           ddudydy = 0.
           ddudydz = 0.
           ddudzdz = 0.
           do i = 1 , ib_stencil

              imp(1,1,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,1,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,1,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));

              imp(1,2,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,2,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,2,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));

              imp(1,3,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,3,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,3,2) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));

              ue   = ue   + phile(i,1)*imp(2,2,2);
              dudx = dudx + phile(i,1)*(imp(3,2,2)-imp(1,2,2))/(2.*dx)
              dudy = dudy + phile(i,1)*(imp(2,3,2)-imp(2,1,2))/(2.*dy)
              ddudxdx = ddudxdx + phile(i,1)*(imp(3,2,2)-2.*imp(2,2,2)+imp(1,2,2))/(dx*dx)
              ddudydy = ddudydy + phile(i,1)*(imp(2,3,2)-2.*imp(2,2,2)+imp(2,1,2))/(dy*dy)
              ddudxdy = ddudxdy + phile(i,1)*(imp(3,3,2)-imp(3,1,2)-imp(1,3,2)+imp(1,1,2))/(4.*dx*dy)
#if NDIM == MDIM

!              imp(1,1,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(2,1,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1)-1);
!              imp(3,1,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);

              imp(1,2,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(2,2,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(3,2,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);

!              imp(1,3,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(2,3,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1)-1);
!              imp(3,3,1) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);
              
              !------------------
!              imp(1,1,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(2,1,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1)+1);
!              imp(3,1,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);

              imp(1,2,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(2,2,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(3,2,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);

!              imp(1,3,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(2,3,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1)+1);
!              imp(3,3,3) = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);

              dudz = dudz + phile(i,1)*(imp(2,2,3)-imp(2,2,1))/(2.*dz)
              ddudzdz = ddudzdz + phile(i,1)*(imp(2,2,3)-2.*imp(2,2,2)+imp(2,2,1))/(dz*dz)
              ddudxdz = ddudxdz + phile(i,1)*(imp(3,2,3)-imp(3,2,1)-imp(1,2,3)+imp(1,2,1))/(4.*dx*dz)
              ddudydz = ddudydz + phile(i,1)*(imp(2,3,3)-imp(2,3,1)-imp(2,1,3)+imp(2,1,1))/(4.*dy*dz)
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facexData,FACEX)

           case(2)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,faceyData,FACEY)

           ! V velocity derivatives on external point:
           ve = 0.; dvdx = 0.; dvdy =0.;
           dvdz = 0.
           ddvdxdx = 0.
           ddvdxdy = 0.
           ddvdxdz = 0.
           ddvdydy = 0.
           ddvdydz = 0.
           ddvdzdz = 0.

           do i = 1 , ib_stencil

              imp(1,1,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,1,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,1,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));

              imp(1,2,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,2,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,2,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));

              imp(1,3,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,3,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,3,2) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));

              ve   = ve   + phile(i,1)*imp(2,2,2);
              dvdx = dvdx + phile(i,1)*(imp(3,2,2)-imp(1,2,2))/(2.*dx)
              dvdy = dvdy + phile(i,1)*(imp(2,3,2)-imp(2,1,2))/(2.*dy)
              ddvdxdx = ddvdxdx + phile(i,1)*(imp(3,2,2)-2.*imp(2,2,2)+imp(1,2,2))/(dx*dx)
              ddvdydy = ddvdydy + phile(i,1)*(imp(2,3,2)-2.*imp(2,2,2)+imp(2,1,2))/(dy*dy)
              ddvdxdy = ddvdxdy + phile(i,1)*(imp(3,3,2)-imp(3,1,2)-imp(1,3,2)+imp(1,1,2))/(4.*dx*dy)

#if NDIM == MDIM

!              imp(1,1,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(2,1,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1)-1);
!              imp(3,1,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);

              imp(1,2,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(2,2,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(3,2,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);

!              imp(1,3,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(2,3,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1)-1);
!              imp(3,3,1) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)-1);
              
              !------------------
!              imp(1,1,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(2,1,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1)+1);
!              imp(3,1,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)-1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);

              imp(1,2,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(2,2,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(3,2,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);

!              imp(1,3,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(2,3,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1)+1);
!              imp(3,3,3) = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
!                                                   ielem(i,JAXIS,presflag+1)+1, &
!                                                   ielem(i,KAXIS,presflag+1)+1);

              dvdz = dvdz + phile(i,1)*(imp(2,2,3)-imp(2,2,1))/(2.*dz)
              ddvdzdz = ddvdzdz + phile(i,1)*(imp(2,2,3)-2.*imp(2,2,2)+imp(2,2,1))/(dz*dz)
              ddvdxdz = ddvdxdz + phile(i,1)*(imp(3,2,3)-imp(3,2,1)-imp(1,2,3)+imp(1,2,1))/(4.*dx*dz)
              ddvdydz = ddvdydz + phile(i,1)*(imp(2,3,3)-imp(2,3,1)-imp(2,1,3)+imp(2,1,1))/(4.*dy*dz)
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


           case(3) 

#if NDIM == MDIM
           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facezData,FACEZ)

           ! W velocity derivatives on external point:
           we   = 0.; dwdx = 0.; dwdy =0.;
           dwdz = 0.
           ddwdxdx = 0.
           ddwdxdy = 0.
           ddwdxdz = 0.
           ddwdydy = 0.
           ddwdydz = 0.
           ddwdzdz = 0.
           do i = 1 , ib_stencil

              imp(2,1,1) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1)-1);

              imp(1,2,1) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(2,2,1) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);
              imp(3,2,1) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)-1);

              imp(2,3,1) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1)-1);
              
              !------------------
              imp(2,1,3) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1)+1);

              imp(1,2,3) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(2,2,3) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);
              imp(3,2,3) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1)+1);

              imp(2,3,3) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1)+1);
              !-----------------------
              imp(1,1,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,1,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,1,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1)-1, &
                                                   ielem(i,KAXIS,presflag+1));

              imp(1,2,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,2,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,2,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1), &
                                                   ielem(i,KAXIS,presflag+1));

              imp(1,3,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)-1, &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(2,3,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));
              imp(3,3,2) = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1)+1, &
                                                   ielem(i,JAXIS,presflag+1)+1, &
                                                   ielem(i,KAXIS,presflag+1));

              we   = we   + phile(i,1)*imp(2,2,2);
              dwdx = dwdx + phile(i,1)*(imp(3,2,2)-imp(1,2,2))/(2.*dx)
              dwdy = dwdy + phile(i,1)*(imp(2,3,2)-imp(2,1,2))/(2.*dy)
              dwdz = dwdz + phile(i,1)*(imp(2,2,3)-imp(2,2,1))/(2.*dz)
              ddwdxdx = ddwdxdx + phile(i,1)*(imp(3,2,2)-2.*imp(2,2,2)+imp(1,2,2))/(dx*dx)
              ddwdydy = ddwdydy + phile(i,1)*(imp(2,3,2)-2.*imp(2,2,2)+imp(2,1,2))/(dy*dy)
              ddwdxdy = ddwdxdy + phile(i,1)*(imp(3,3,2)-imp(3,1,2)-imp(1,3,2)+imp(1,1,2))/(4.*dx*dy)
              ddwdzdz = ddwdzdz + phile(i,1)*(imp(2,2,3)-2.*imp(2,2,2)+imp(2,2,1))/(dz*dz)
              ddwdxdz = ddwdxdz + phile(i,1)*(imp(3,2,3)-imp(3,2,1)-imp(1,2,3)+imp(1,2,1))/(4.*dx*dz)
              ddwdydz = ddwdydz + phile(i,1)*(imp(2,3,3)-imp(2,3,1)-imp(2,1,3)+imp(2,1,1))/(4.*dy*dz)

           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

           end select

#endif
           
        end if

     enddo

  enddo

!============================
! Project external and marker Velocity on the plane:
  ven = nxp*ue + nyp*ve + nzp*we

  ! ves = ve - (ve*n) n
  ues  = ue - ven*nxp
  ves  = ve - ven*nyp
  wes  = we - ven*nzp

  ! vps = vp - (vp*n) n
  vpn  = nxp*ubd + nyp*vbd + nzp*wbd
  ups  = ubd - vpn*nxp
  vps  = vbd - vpn*nyp
  wps  = wbd - vpn*nzp

  ! Linear Part:
  normt = 1./h
  dun = (ues-ups)*normt
  dvn = (ves-vps)*normt
  dwn = (wes-wps)*normt
 
  ! First versor in the local tangent velocity direction:
  if( (abs(dun)+abs(dvn)+abs(dwn)) .lt. 1.0e-14) then ! case zero velocity difference
    tp = 0.
  else
   normt = 1./sqrt(dun**2. + dvn**2. + dwn**2.)
   tp(1) = dun*normt
   tp(2) = dvn*normt
   tp(3) = dwn*normt
  endif


  exx = dudx
  exy = 0.5*(dudy + dvdx)
  exz = 0.5*(dudz + dwdx)
  eyy = dvdy
  eyz = 0.5*(dvdz + dwdy)
  ezz = dwdz

  ddvel(1,1,1) = ddudxdx
  ddvel(2,1,1) = ddvdxdx
  ddvel(3,1,1) = ddwdxdx
  ddvel(1,2,1) = ddudxdy
  ddvel(2,2,1) = ddvdxdy
  ddvel(3,2,1) = ddwdxdy
  ddvel(1,3,1) = ddudxdz
  ddvel(2,3,1) = ddvdxdz
  ddvel(3,3,1) = ddwdxdz

  ddvel(1,1,2) = ddudxdy
  ddvel(2,1,2) = ddvdxdy
  ddvel(3,1,2) = ddwdxdy
  ddvel(1,2,2) = ddudydy
  ddvel(2,2,2) = ddvdydy
  ddvel(3,2,2) = ddwdydy
  ddvel(1,3,2) = ddudydz
  ddvel(2,3,2) = ddvdydz
  ddvel(3,3,2) = ddwdydz

  ddvel(1,1,3) = ddudxdz
  ddvel(2,1,3) = ddvdxdz
  ddvel(3,1,3) = ddwdxdz
  ddvel(1,2,3) = ddudydz
  ddvel(2,2,3) = ddvdydz
  ddvel(3,2,3) = ddwdydz
  ddvel(1,3,3) = ddudzdz
  ddvel(2,3,3) = ddvdzdz
  ddvel(3,3,3) = ddwdzdz

  ddp(1,1) = ddpdxdx
  ddp(2,1) = ddpdxdy
  ddp(3,1) = ddpdxdz
  ddp(1,2) = ddpdxdy
  ddp(2,2) = ddpdydy
  ddp(3,2) = ddpdydz
  ddp(1,3) = ddpdxdz
  ddp(2,3) = ddpdydz
  ddp(3,3) = ddpdzdz


  dpdxn = dpdx*nxp + dpdy*nyp + dpdz*nzp
  zL2    = zpres - dpdxn*h

  ! Fvisc = Tau * n
  fvx_l = nu*2.*(exx*nxp+exy*nyp+exz*nzp)
  fvy_l = nu*2.*(exy*nxp+eyy*nyp+eyz*nzp)
  fvz_l = nu*2.*(exz*nxp+eyz*nyp+ezz*nzp)

  fv_c = 0.0
  do kk = 1, MDIM
    do jj = 1, MDIM
      do ii = 1, MDIM
        fv_c = fv_c + tp(ii)*np(jj)*np(kk)*ddvel(ii,jj,kk)
      enddo
    enddo
  enddo  

  fv_c  =  -fv_c*h*nu
  fvx_c = fv_c*tp(1)
  fvy_c = fv_c*tp(2)
  fvz_c = fv_c*tp(3)

  fv_cc = 0.0
  do jj = 1, MDIM
    do ii = 1, MDIM
      fv_cc = fv_cc +  tp(ii)*np(jj)*ddp(ii,jj)
    enddo
  enddo  
  fv_cc = fv_cc*h*h*0.5   ! Notice: There is no nu here

  fvx_cc = fv_cc*tp(1)
  fvy_cc = fv_cc*tp(2)
  fvz_cc = fv_cc*tp(3)

!  w_l(1) = (fvy_l*nxp - fvx_l*nyp)/nu
!  w_l(2) = (fvz_l*nyp - fvy_l*nzp)/nu
!  w_l(3) = (fvx_l*nzp - fvz_l*nxp)/nu

  w_l(1) = dvdx - dudy ! w_z
  w_l(2) = dwdy - dvdz ! w_x
  w_l(3) = dudz - dwdx ! w_y

  w_c(1) = (fvy_c*nxp - fvx_c*nyp)/nu
  w_c(2) = (fvz_c*nyp - fvy_c*nzp)/nu
  w_c(3) = (fvx_c*nzp - fvz_c*nxp)/nu

  w_cc(1) = (fvy_cc*nxp - fvx_cc*nyp)/nu
  w_cc(2) = (fvz_cc*nyp - fvy_cc*nzp)/nu
  w_cc(3) = (fvx_cc*nzp - fvz_cc*nxp)/nu

  zv = w_l + w_c + w_cc

  particleData(FXVI_PART_PROP) = fvx_l + fvx_c + fvx_cc
  particleData(FYVI_PART_PROP) = fvy_l + fvy_c + fvy_cc
  particleData(FZVI_PART_PROP) = fvz_l + fvz_c + fvz_cc

  particleData(PRES_PART_PROP) = zL2
  particleData(PEX0_PART_PROP) = particleData(PEXT_PART_PROP)
  particleData(PEXT_PART_PROP) = zpres

#ifdef TEST_COMPARE
  write(9800+gr_meshMe,'(10f20.12)') particleData(GLOB_PART_PROP), zL2, zL, zpres 
  write(9900+gr_meshMe,'(30f20.12)') particleData(GLOB_PART_PROP), zv(1), w_l(1), w_c(1), w_cc(1) 
!  write(9700+gr_meshMe,'(30f20.12)') particleData(GLOB_PART_PROP), ddpdxdx,ddpdxdy,ddpdydy
#endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else /* NORMAL_GRAD_CORR */

  zpres = 0.
  zv(1:MDIM) = 0.
  zL= 0.
  nuwx = 0.
  nuwy = 0.
  nuwz = 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
#ifdef TANGENT_WITH_VORTICITY
     nkij = 1 + (1-presflag)*(2*NDIM-MDIM-1)
#else
     nkij = 1 + (1-presflag)*(NDIM-1) ! in 2D d/dx,d/dy; in 3D d/dx,d/dy,d/dz
#endif

     do gridind = 1,nkij

#ifdef TANGENT_WITH_VORTICITY
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdnu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlnu(idim,gridind))*del(idim)
        enddo
#else
        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdu(:,gridind)
        ! Auxiliary deltas
        do idim = 1,NDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlu(idim,gridind))*del(idim)
        enddo
#endif

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe,np,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe,xyz_stencil,del,derivflag,phile)
    
        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

#ifndef TANGENT_WITH_VORTICITY
           dpdx = 0.; dpdy =0.;
           dpdz = 0.
#endif
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                      ielem(i,JAXIS,presflag+1), &
                                      ielem(i,KAXIS,presflag+1));  

              zpres = zpres + phile(i,1)*p_i 

#ifndef TANGENT_WITH_VORTICITY
              dpdx = dpdx + phile(i,2)*p_i;
              dpdy = dpdy + phile(i,3)*p_i;
#if NDIM == MDIM
              dpdz = dpdz + phile(i,4)*p_i;
#endif             
#endif
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           ! Get Pressure approximation at surface marker: Acceleration +
           ! gravity effects.
           dpdn = -(     ubdd*nxp +      vbdd*nyp +      wbdd*nzp) + & ! -rho*Du/Dt * n 
                   (ins_gravX*nxp + ins_gravY*nyp + ins_gravZ*nzp);    ! +rho*    g * n
           zL = zpres - dpdn*h;

        else                                         ! Tangent stress

#ifdef TANGENT_WITH_VORTICITY
 
           select case (gridind)
           case(1) 
           ! Only wz:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortz(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwz = nu*zv(gridind)

           case(2)
           ! wx:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortx(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwx = nu*zv(gridind)          

           case(3)
           ! wy:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vorty(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwy = nu*zv(gridind)

           end select

#else

           select case(gridind)
           case(1)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facexData,FACEX)

           ! U velocity derivatives on external point:
           ue = 0.; dudx = 0.; dudy =0.;
           dudz = 0.
           do i = 1 , ib_stencil

              ui = facexData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              ue   = ue   + phile(i,1)*ui;
              dudx = dudx + phile(i,2)*ui;
              dudy = dudy + phile(i,3)*ui;
#if NDIM == MDIM
              dudz = dudz + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facexData,FACEX)

           case(2)

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,faceyData,FACEY)

           ! V velocity derivatives on external point:
           ve = 0.; dvdx = 0.; dvdy =0.;
           dvdz = 0.
           do i = 1 , ib_stencil

              ui = faceyData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              ve   = ve   + phile(i,1)*ui; 
              dvdx = dvdx + phile(i,2)*ui;
              dvdy = dvdy + phile(i,3)*ui;
#if NDIM == MDIM
              dvdz = dvdz + phile(i,4)*ui;
#endif             
           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


           case(3) 

#if NDIM == MDIM
           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,facezData,FACEZ)

           ! W velocity derivatives on external point:
           we   = 0.; dwdx = 0.; dwdy =0.;
           dwdz = 0.
           do i = 1 , ib_stencil

              ui = facezData(VELC_FACE_VAR,ielem(i,IAXIS,presflag+1), &
                                           ielem(i,JAXIS,presflag+1), &
                                           ielem(i,KAXIS,presflag+1));

              we   = we   + phile(i,1)*ui;
              dwdx = dwdx + phile(i,2)*ui;
              dwdy = dwdy + phile(i,3)*ui;
              dWdz = dwdz + phile(i,4)*ui;

           enddo

           ! Release pointers:
           call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

           end select

#endif
           
        end if

     enddo

  enddo

  ! Assign pressure and viscous forces to particleData:
#ifdef INS_CONSTDENS
  ! Case constant density:
  particleData(PRES_PART_PROP) = zL !+ ins_gravX*(xp-gr_imin) + &
                                    !  ins_gravY*(yp-gr_jmin) + &
                                    !  ins_gravZ*(zp-gr_kmin)
  particleData(PEX0_PART_PROP) = particleData(PEXT_PART_PROP)
  particleData(PEXT_PART_PROP) = zpres !+ ins_gravX*(xbe(IAXIS)-gr_imin) + &
                                       !  ins_gravY*(xbe(JAXIS)-gr_jmin) + &
                                       !  ins_gravZ*(xbe(KAXIS)-gr_kmin)


#ifdef TANGENT_WITH_VORTICITY 

  ! Fvisc = -nu (N x W)
  particleData(FXVI_PART_PROP) = (nzp*nuwy-nyp*nuwz)
  particleData(FYVI_PART_PROP) = (nxp*nuwz-nzp*nuwx)
  particleData(FZVI_PART_PROP) = (nyp*nuwx-nxp*nuwy)

#else

#ifdef USE_CF

  ! Project external and marker Velocity on the plane:
  ven = nxp*ue + nyp*ve + nzp*we

  ! ves = ve - (ve*n) n
  ues  = ue - ven*nxp
  ves  = ve - ven*nyp
  wes  = we - ven*nzp

  ! vps = vp - (vp*n) n
  vpn  = nxp*ubd + nyp*vbd + nzp*wbd
  ups  = ubd - vpn*nxp
  vps  = vbd - vpn*nyp
  wps  = wbd - vpn*nzp

  ! Linear Part:
  normt = 1./h
  dun = (ues-ups)*normt
  dvn = (ves-vps)*normt
  dwn = (wes-wps)*normt

  ! Correction from diffusion equation (A. Posa):
  ! First versor in the local tangent velocity direction:
  if( (abs(dun)+abs(dvn)+abs(dwn)) .lt. 1.0e-14) then ! case zero velocity difference
  tx = 0.
  ty = 0.
  tz = 0.
  dpdxt = 0.
  else
  normt = 1./sqrt(dun**2. + dvn**2. + dwn**2.)
  tx = dun*normt
  ty = dvn*normt
  tz = dwn*normt

  ! Compute dpdxt from pressure gradients directly, use only hydrodynamic dpdxt:
  dpdxt = (dpdx-ins_gravX)*tx + (dpdy-ins_gravY)*ty + (dpdz-ins_gravZ)*tz

#ifdef TWO_POINTSP
  ! Compute dpdxt from other pressure values:
  ! External Point Position:
  normt = h/2.
  xbe2(IAXIS) = xp + nxp*h + tx*normt
  xbe2(JAXIS) = yp + nyp*h + ty*normt
#if NDIM == 3
  xbe2(KAXIS) = zp + nzp*h + tz*normt
#else
  xbe2(KAXIS) = 0.
#endif
  zpres2 = 0.
  presflag = CONSTANT_ONE
  gridind  = 1
  ! Define Grids in case of pressure:
  gridfl(:) = grdip(:)
  ! Auxiliary deltas
  do idim = 1,MDIM
     delaux(idim) = (dlip(idim))*del(idim)
  enddo

  !! Point 2:
  ! Obtain Stencil for External Point:
  call ib_stencils(xbe2,np,gridfl,del,coord,bsize,   &
                   ielem(:,:,presflag+1),hl,COMPUTE_FORCES)
  ! Compute shape functions
  ! Positions of points on the stencil:
  xyz_stencil(1:ib_stencil,1:MDIM) = 0.
  do idim = 1,NDIM
     xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim)
  enddo
  ! Get interpolation functions:
  call ib_getInterpFunc(xbe2,xyz_stencil,del,0,phile)
  ! Point to cell centered Variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  ! Value of the function in xbe:
  do i = 1 , ib_stencil
      p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                              ielem(i,JAXIS,presflag+1), &
                              ielem(i,KAXIS,presflag+1));
      zpres2 = zpres2 + phile(i,1)*p_i
  enddo
  ! Release Pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  ! Substract hydrostatic pressure to total pressure in point 2:
  zpres2 = zpres2 - ins_gravX*(xbe2(IAXIS)-gr_imin) - &
                    ins_gravY*(xbe2(JAXIS)-gr_jmin) - &
                    ins_gravZ*(xbe2(KAXIS)-gr_kmin)

  !! Point 3:
  xbe3(IAXIS) = xp + nxp*h - tx*normt
  xbe3(JAXIS) = yp + nyp*h - ty*normt
#if NDIM == 3
  xbe3(KAXIS) = zp + nzp*h - tz*normt
#else
  xbe3(KAXIS) = 0.
#endif
  zpres3 = 0.
  ! Obtain Stencil for External Point:
  call ib_stencils(xbe3,np,gridfl,del,coord,bsize,   &
                   ielem(:,:,presflag+1),hl,COMPUTE_FORCES)
  ! Compute shape functions
  ! Positions of points on the stencil:
  xyz_stencil(1:ib_stencil,1:MDIM) = 0.
  do idim = 1,NDIM
     xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim)
  enddo
  ! Get interpolation functions:
  call ib_getInterpFunc(xbe3,xyz_stencil,del,0,phile)
  ! Point to cell centered Variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  ! Value of the function in xbe:
  do i = 1 , ib_stencil
      p_i = solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                              ielem(i,JAXIS,presflag+1), &
                              ielem(i,KAXIS,presflag+1));
      zpres3 = zpres3 + phile(i,1)*p_i
  enddo
  ! Release Pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  ! Substract hydrostatic pressure to total pressure in point 3:
  zpres3 = zpres3 - ins_gravX*(xbe3(IAXIS)-gr_imin) - &
                    ins_gravY*(xbe3(JAXIS)-gr_jmin) - &
                    ins_gravZ*(xbe3(KAXIS)-gr_kmin)

  ! Finally dpdxt:
  !dpdxt = (zpres2-zpres)/normt ! This 1st order approx would save us from
                                ! computing zpres3
  dpdxt = (zpres2-zpres3)/(2.*normt) 

#endif /* TWO POINT flag */

  endif ! case zero velocity difference

  ! Diffusion correction - 1/rho*dp/dxt*h/(2nu) * t:
  normt = h/(2.*nu)

  ! Using Strain tensor on the external point:
  ! Strain velocities tensor:
  exx = dudx
  exy = 0.5*(dudy + dvdx)
  exz = 0.5*(dudz + dwdx)
  eyy = dvdy
  eyz = 0.5*(dvdz + dwdy)
  ezz = dwdz

  ! Fvisc = Tau * n
  dun = 2.*(exx*nxp+exy*nyp+exz*nzp) - (dpdxt*tx)*normt
  dvn = 2.*(exy*nxp+eyy*nyp+eyz*nzp) - (dpdxt*ty)*normt
  dwn = 2.*(exz*nxp+eyz*nyp+ezz*nzp) - (dpdxt*tz)*normt

  ! Using Linear du/dxn:
!  dun = dun - (dpdxt*tx)*normt
!  dvn = dvn - (dpdxt*ty)*normt
!  dwn = dwn - (dpdxt*tz)*normt
 
  ! With quadratic correction:
  !dune = dudx*nxp + dudy*nyp + dudz*nzp
  !dvne = dvdx*nxp + dvdy*nyp + dvdz*nzp
  !dwne = dwdx*nxp + dwdy*nyp + dwdz*nzp 
  !dun = 2.*(ues-ups)/h - dune
  !dvn = 2.*(ves-vps)/h - dvne
  !dwn = 2.*(wes-wps)/h - dwne

  ! Fvisc = nu dv/dn - 1/rho*dp/dxt*h/2 * t, here rho=1:
  particleData(FXVI_PART_PROP) = nu*dun 
  particleData(FYVI_PART_PROP) = nu*dvn 
  particleData(FZVI_PART_PROP) = nu*dwn 
 
  ! Vorticity:
  ! In Z dir: wz = dv/dx - du/dy
  zv(1) = dvn*nxp - dun*nyp
  ! In X dir: wx = dw/dy - dv/dz
  zv(2) = dwn*nyp - dvn*nzp
  ! In Y dir: wy = du/dz - dw/dx
  zv(3) = dun*nzp - dwn*nxp


#else /* USE_CF */


  ! Strain velocities tensor:
  exx = dudx
  exy = 0.5*(dudy + dvdx)
  exz = 0.5*(dudz + dwdx)
  eyy = dvdy 
  eyz = 0.5*(dvdz + dwdy)
  ezz = dwdz 

  ! Fvisc = Tau * n
  particleData(FXVI_PART_PROP) = 2.*nu*(exx*nxp+exy*nyp+exz*nzp)
  particleData(FYVI_PART_PROP) = 2.*nu*(exy*nxp+eyy*nyp+eyz*nzp)
  particleData(FZVI_PART_PROP) = 2.*nu*(exz*nxp+eyz*nyp+ezz*nzp)

  ! Vorticity:
  ! In Z dir: wz = dv/dx - du/dy
  zv(1) = dvdx - dudy
  ! In X dir: wx = dw/dy - dv/dz
  zv(2) = dwdy - dvdz
  ! In Y dir: wy = du/dz - dw/dx
  zv(3) = dudz - dwdx


#endif
#endif

#else /* not INS_CONSTDENS */

  call Driver_abortFlash("Variable density particle pressure and tangent stress not defined.")

#endif

#endif /* NORMAL_GRAD_CORR */

  ! Vorticity components:
  particleData(VORZ_PART_PROP) = zv(1)
  particleData(VORX_PART_PROP) = zv(2)
  particleData(VORY_PART_PROP) = zv(3)

  return

end subroutine ib_distributedForces


