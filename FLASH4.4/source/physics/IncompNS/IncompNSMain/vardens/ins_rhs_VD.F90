
      SUBROUTINE ins_rhs2d_VD(uni,vni,ru1,ix1,ix2,jy1,jy2,dx,dy,ru,rv, &
                           visc,rho1x,rho2x,rho1y,rho2y,gravX,gravY)

  !***************************************************************
  ! This subroutine computes the discretization of the RHS of the 
  ! Helmholtz equation on a staggered uniform grid.
  !
  ! Input:  uni,vni     = velocity at timestep n
  !         ru1         = molecular viscosity !- kpd - Inverse Reynolds No
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         dx,dy       = grid spacing in x and y directions
  !
  ! Output: ru,rv    = u and v momentum for Helmholtz RHS
  !**************************************************************

      use Driver_interface, ONLY : Driver_abortFlash

      use Driver_data,      ONLY : dr_nstep

      use RuntimeParameters_interface, ONLY : RuntimeParameters_get

      use IncompNS_data, ONLY : ins_iConvU

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
      REAL, INTENT(IN):: ru1, dx, dy, gravX,gravY
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, visc, rho1x, rho2x, rho1y, rho2y
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv

      INTEGER:: i, j
      REAL:: dx1, dy1, Mdens
      ! x-component variables
      REAL:: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
             uyplus, uyminus
      REAL:: dudxp, dudxm, dudyp, dudym, dvdxp, dvdxm
      REAL:: tvjp, tvjm
      REAL:: txxp, txxm, tyyp, tyym
      REAL:: txyp, txym
      ! new y-component variables
      REAL:: vyplus, vyminus
      REAL:: dvdyp, dvdym
      REAL:: tvip, tvim
      INTEGER, parameter :: kz1 = 1

      !KPD - 1st or 3rd order upwind
      REAL :: uu,uxp,uxm,uyp,uym,vv,vxp,vxm,vyp,vym,uy,vx,ududx,udvdx,vdudy,vdvdy
      REAL :: uxpp,uxmm,uypp,uymm,vxpp,vxmm,vypp,vymm
      INTEGER :: iConvU
      REAL :: rConvU

      !KPD - Choose Convection Scheme 1=1st order upwind, 2=2nd order central, 3=3rd order upwind
      !iConvU = 4
      iConvU = ins_iConvU
      !print*,"iConvU VALUE: ",iConvU

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy

      !++++++++++  U-COMPONENT  ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2+1

             !=============================================================
             !KPD - 1st Order Upwind... ===================================
             uu   = uni(i,j,kz1)

             uxp  = uni(i+1,j,kz1)
             uxm  = uni(i-1,j,kz1)
             uxpp = uni(i+2,j,kz1)
             uxmm = uni(i-2,j,kz1)

             uyp  = uni(i,j+1,kz1)
             uym  = uni(i,j-1,kz1)
             uypp = uni(i,j+2,kz1)
             uymm = uni(i,j-2,kz1)

             uy   = 0.25*( uni(i,j,kz1) + uni(i+1,j,kz1) + uni(i,j-1,kz1) + uni(i+1,j-1,kz1) ) 

             vv   = vni(i,j,kz1)

             vxp  = vni(i+1,j,kz1)
             vxm  = vni(i-1,j,kz1)
             vxpp = vni(i+2,j,kz1)
             vxmm = vni(i-2,j,kz1)

             vyp  = vni(i,j+1,kz1)
             vym  = vni(i,j-1,kz1)
             vypp = vni(i,j+2,kz1)
             vymm = vni(i,j-2,kz1)

             vx  = 0.25*( vni(i,j,kz1) + vni(i-1,j,kz1) + vni(i,j+1,kz1) + vni(i-1,j+1,kz1) )

             !=============================================================
             ! u.grad(u) = uj*dui/dxj = [ u*du/dx + v*du/dy ] (i) + [ u*dv/dx + v*dv/dy ] (j)
             !=============================================================
             if (uu .gt. 0) then
                if (iConvU .eq. 1) then
                   ududx = uu*(uu-uxm)*dx1   
                elseif (iConvU .eq. 2) then
                   ududx = uu*(3.*uu - 4.*uxm + 1.*uxmm)/2.*dx1   
                elseif (iConvU .eq. 3) then
                   ududx = uu*(2.*uxp + 3.*uu - 6.*uxm + 1.*uxmm)/6.*dx1   
                end if
             else
                if (iConvU .eq. 1) then
                  !ududx = uu*(uu-uxp)*dx1   
                   ududx = uu*(uxp-uu)*dx1   
                elseif (iConvU .eq. 2) then
                   ududx = uu*(-3.*uu + 4.*uxp - 1.*uxpp)/2.*dx1   
                elseif (iConvU .eq. 3) then
                  !ududx = uu*( 2.*uxm + 3.*uu - 6.*uxp + 1.*uxpp)/6.*dx1   
                   ududx = uu*(-2.*uxm - 3.*uu + 6.*uxp - 1.*uxpp)/6.*dx1   
                end if
             end if
             if (vx .gt. 0) then
                if (iConvU .eq. 1) then
                   vdudy = vx*(uu-uym)*dy1 
                elseif (iConvU .eq. 2) then
                   vdudy = vx*(3.*uu - 4.*uym + 1.*uymm)/2.*dy1 
                elseif (iConvU .eq. 3) then
                   vdudy = vx*(2.*uyp +3.*uu - 6.*uym + 1.*uymm)/6.*dy1 
                end if
             else
                if (iConvU .eq. 1) then
                  !vdudy = vx*(uu-uyp)*dy1
                   vdudy = vx*(uyp-uu)*dy1
                elseif (iConvU .eq. 2) then
                   vdudy = vx*(-3.*uu + 4.*uyp - 1.*uypp)/2.*dy1 
                elseif (iConvU .eq. 3) then
                  !vdudy = vx*( 2.*uym +3.*uu - 6.*uyp + 1.*uypp)/6.*dy1 
                   vdudy = vx*(-2.*uym -3.*uu + 6.*uyp - 1.*uypp)/6.*dy1 
                end if
             end if
             !=============================================================
             !=============================================================

             ! get velocities at 1/2 locations
             uxplus = (uni(i+1,j,kz1) + uni(i,j,kz1))*0.5
             uxminus = (uni(i,j,kz1) + uni(i-1,j,kz1))*0.5

             vxplus = (vni(i,j+1,kz1) + vni(i-1,j+1,kz1))*0.5
             vxminus = (vni(i,j,kz1) + vni(i-1,j,kz1))*0.5

             uyplus = (uni(i,j+1,kz1) + uni(i,j,kz1))*0.5
             uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             ! get derivatives at 1/2 locations
             dudxp = (uni(i+1,j,kz1) - uni(i,j,kz1))*dx1
             dudxm = (uni(i,j,kz1) - uni(i-1,j,kz1))*dx1
             dudyp = (uni(i,j+1,kz1) - uni(i,j,kz1))*dy1
             dudym = (uni(i,j,kz1) - uni(i,j-1,kz1))*dy1
             dvdxp = (vni(i,j+1,kz1) - vni(i-1,j+1,kz1))*dx1    !- kpd - Not Used?
             dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1))*dx1        !- kpd - Not Used?

             ! flux of normal total stresses
             !------------------------------

             !- kpd - Constant Viscosity Implementation
             !txxp = ru1*dudxp                !- kpd - should use VISC(i  ,j   )
             !txxm = ru1*dudxm                !- kpd - should use VISC(i-1,j   )
             !tyyp = ru1*dudyp                !- kpd - should use VISC(i  ,j+.5) 
             !tyym = ru1*dudym                !- kpd - should use VISC(i  ,j-.5)

             if (visc(i,j,kz1) .lt. 0.0 .OR. visc(i,j,kz1) .gt. 1.0) then
                print*,"ERROR1: Bad Viscosity Value in the RHS",i,j,visc(i,j,kz1)
                call Driver_abortFlash('ERROR: Bad RHS Viscosity')
             end if

             !- kpd - Variable Viscosity Implementation (ru1 is 1/Re)
             txxp = ru1*visc(i  ,j,kz1)*dudxp                   
             txxm = ru1*visc(i-1,j,kz1)*dudxm                      
             tyyp = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i-1,j+1,kz1)+visc(i,j+1,kz1))*dudyp    
             tyym = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i-1,j-1,kz1)+visc(i,j-1,kz1))*dudym   

             Mdens = ( rho1x(i,j,kz1) + rho2x(i,j,kz1) )  ! Mixture inverse density

             if (Mdens .lt. 1.0 .OR. Mdens .gt. 5000.0) then
                print*,"ERROR: Bad Density Value in the RHS",i,j,Mdens
                call Driver_abortFlash('ERROR: Bad RHS Viscosity')
             end if


             if (iConvU .eq. 0) then
                rConvU = - (uxplus*uxplus - uxminus*uxminus)*dx1   &! advection term
                         - (vxplus*uyplus - vxminus*uyminus)*dy1
             elseif (iConvU .eq. 1 .OR. iConvU .eq. 3 .OR. iConvU .eq. 2) then
                rConvU = - ududx - vdudy 
             end if

             ! calculate RHS for u-momentum
             ru(i,j,kz1) =                                          &
                         !- (uxplus*uxplus - uxminus*uxminus)*dx1   &! advection term
                         !- (vxplus*uyplus - vxminus*uyminus)*dy1   &
                         !- ududx                                   &! advection term
                         !- vdudy                                   & 
                          rConvU                                    &
                          + Mdens*(txxp - txxm)*dx1                 & ! diffusion - normal terms 
                          + Mdens*(tyyp - tyym)*dy1                 &
                          - gravX     
             !!ru(i,j,kz1) = 0.0 
          enddo
       enddo


    !++++++++++  V-COMPONENT  ++++++++++

       do j = jy1,jy2+1
          do i = ix1,ix2

             !=============================================================
             !KPD - 1st Order Upwind... ===================================
             uu   = uni(i,j,kz1)

             uxp  = uni(i+1,j,kz1)
             uxm  = uni(i-1,j,kz1)
             uxpp = uni(i+2,j,kz1)
             uxmm = uni(i-2,j,kz1)

             uyp  = uni(i,j+1,kz1)
             uym  = uni(i,j-1,kz1)
             uypp = uni(i,j+2,kz1)
             uymm = uni(i,j-2,kz1)

             uy   = 0.25*( uni(i,j,kz1) + uni(i+1,j,kz1) + uni(i,j-1,kz1) + uni(i+1,j-1,kz1) ) 

             vv   = vni(i,j,kz1)

             vxp  = vni(i+1,j,kz1)
             vxm  = vni(i-1,j,kz1)
             vxpp = vni(i+2,j,kz1)
             vxmm = vni(i-2,j,kz1)

             vyp  = vni(i,j+1,kz1)
             vym  = vni(i,j-1,kz1)
             vypp = vni(i,j+2,kz1)
             vymm = vni(i,j-2,kz1)

             vx  = 0.25*( vni(i,j,kz1) + vni(i-1,j,kz1) + vni(i,j+1,kz1) + vni(i-1,j+1,kz1) )


             !=============================================================
             ! u.grad(u) = uj*dui/dxj = [ u*du/dx + v*du/dy ] (i) + [ u*dv/dx + v*dv/dy ] (j)
             !=============================================================
             if (uy .gt. 0) then
                if (iConvU .eq. 1) then
                   udvdx = uy*(vv-vxm)*dx1 
                elseif (iConvU .eq. 2) then
                   udvdx = uy*(3.*vv - 4.*vxm + 1.*vxmm)/2.*dx1 
                elseif (iConvU .eq. 3) then
                  udvdx = uy*( 2.*vxp + 3.*vv - 6.*vxm + 1.*vxmm)/6.*dx1 
                end if
             else
                if (iConvU .eq. 1) then
                  !udvdx = uy*(vv-vxp)*dx1 
                   udvdx = uy*(vxp-vv)*dx1 
                elseif (iConvU .eq. 2) then
                   udvdx = uy*(-3.*vv + 4.*vxp - 1.*vxpp)/2.*dx1 
                elseif (iConvU .eq. 3) then
                  !udvdx = uy*( 2.*vxm + 3.*vv - 6.*vxp + 1.*vxpp)/6.*dx1 
                   udvdx = uy*(-2.*vxm - 3.*vv + 6.*vxp - 1.*vxpp)/6.*dx1 
                end if
             end if
             if (vv .gt. 0) then
                if (iConvU .eq. 1) then
                   vdvdy = vv*(vv-vym)*dy1 
                elseif (iConvU .eq. 2) then
                   vdvdy = vv*(3.*vv - 4.*vym + 1.*vymm)/2.*dy1 
                elseif (iConvU .eq. 3) then
                   vdvdy = vv*(2.*vyp + 3.*vv - 6.*vym + 1.*vymm)/6.*dy1 
                end if
             else
                if (iConvU .eq. 1) then
                  !vdvdy = vv*(vv-vyp)*dy1 
                   vdvdy = vv*(vyp-vv)*dy1 
                elseif (iConvU .eq. 2) then
                   vdvdy = vv*(-3.*vv + 4.*vyp - 1.*vypp)/2.*dy1 
                elseif (iConvU .eq. 3) then
                  !vdvdy = vv*( 2.*vym + 3.*vv - 6.*vyp + 1.*vypp)/6.*dy1 
                   vdvdy = vv*(-2.*vym - 3.*vv + 6.*vyp - 1.*vypp)/6.*dy1 
                end if
             end if
             !=============================================================
             !=============================================================

             ! get velocities at 1/2 locations
             vxplus = (vni(i+1,j,kz1) + vni(i,j,kz1))*0.5
             vxminus = (vni(i,j,kz1) + vni(i-1,j,kz1))*0.5

             vyplus = (vni(i,j+1,kz1) + vni(i,j,kz1))*0.5
             vyminus = (vni(i,j,kz1) + vni(i,j-1,kz1))*0.5

             uyplus = (uni(i+1,j,kz1) + uni(i+1,j-1,kz1))*0.5
             uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             ! get derivatives at 1/2 locations
             dvdxp = (vni(i+1,j,kz1) - vni(i,j,kz1))*dx1
             dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1))*dx1
             dvdyp = (vni(i,j+1,kz1) - vni(i,j,kz1))*dy1
             dvdym = (vni(i,j,kz1) - vni(i,j-1,kz1))*dy1
             dudyp = (uni(i+1,j,kz1) - uni(i+1,j-1,kz1))*dy1   !- kpd - Not Used?
             dudym = (uni(i,j,kz1) - uni(i,j-1,kz1))*dy1       !- kpd - Not Used?

             ! flux of normal total stresses
             !------------------------------

             !- kpd - Constant Viscosity Implementation
             !txxp = ru1*dvdxp                !- kpd - should use VISC(i+.5,j   )
             !txxm = ru1*dvdxm                !- kpd - should use VISC(i-.5,j   )
             !tyyp = ru1*dvdyp                !- kpd - should use VISC(i   ,j-1 )
             !tyym = ru1*dvdym                !- kpd - should use VISC(i   ,j   )

             !- kpd - Variable Viscosity Implementation (ru1 is 1/Re)
             txxp = ru1*0.25*(visc(i,j,kz1)+visc(i+1,j,kz1)+visc(i,j-1,kz1)+visc(i+1,j-1,kz1))*dvdxp   !- kpd - should use VISC(i+.5,j   )
             txxm = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i,j-1,kz1)+visc(i-1,j-1,kz1))*dvdxm   !- kpd - should use VISC(i-.5,j   )
             tyyp = ru1*visc(i,j,kz1)  *dvdyp                                                          !- kpd - should use VISC(i   ,j-1 )
             tyym = ru1*visc(i,j-1,kz1)*dvdym                                                          !- kpd - should use VISC(i   ,j   )

             Mdens = ( rho1y(i,j,kz1) + rho2y(i,j,kz1) )  ! Mixture inverse density.

             if (Mdens .lt. 1.0 .OR. Mdens .gt. 5000.0) then
                print*,"ERROR: Bad Density Value in the RHS",i,j,Mdens                
                call Driver_abortFlash('ERROR: Bad RHS Viscosity')
             end if

             if (iConvU .eq. 0) then
                rConvU = - (uyplus*vxplus - uyminus*vxminus)*dx1    &! advection term
                         - (vyplus*vyplus - vyminus*vyminus)*dy1 
                         
             elseif (iConvU .eq. 1 .OR. iConvU .eq. 3 .OR. iConvU .eq. 2) then
                rConvU = - udvdx - vdvdy 
             end if
                                                          !         =======
             ! calculate RHS for v-momentum
             rv(i,j,kz1) =                                           &
                         !- (uyplus*vxplus - uyminus*vxminus)*dx1    &! advection term
                         !- (vyplus*vyplus - vyminus*vyminus)*dy1    &
                         !- udvdx                                    &! advection term
                         !- vdvdy                                    & 
                          rConvU                                     &
                          + Mdens* (txxp - txxm)*dx1                 &! diffusion - normal terms
                          + Mdens* (tyyp - tyym)*dy1                 &
                          - gravY                                      ! kpd - gravity term 
             !rv(i,j,kz1) = -grav
          enddo
       enddo


       END SUBROUTINE ins_rhs2d_VD

!*****************************************************************************************************
!*****************************************************************************************************
!*****************************************************************************************************

      SUBROUTINE ins_rhs3d_VD(uni,vni,wni,tv,ru1,      &
                           ix1,ix2,jy1,jy2,kz1,kz2, &
                           dx,dy,dz,ru,rv,rw,visc,  &
                           rho1x,rho2x,rho1y,rho2y, &
                           rho1z,rho2z,gravX, gravY, gravZ)

  !*****************************************************************
  ! This subroutine computes the centered discretization of the RHS 
  ! of the momentum equation (advection + viscous terms) on a 
  ! staggered uniform grid based on the Paramesh grid structure.
  !
  ! Input:  uni,vni,wni = velocity at timestep n
  !         tv          = eddy viscosity
  !         ru1         = molecular viscosity !- kpd - Inverse Reynolds No
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         kz1,kz2     = starting and ending z indices
  !         dx,dy,dz    = grid spacing in x, y, and z directions
  !
  ! Output: ru,rv,rw    = RHS of u, v, and w momentum equations
  !
  ! E. Balaras   July 1999
  ! P. Rabenold  August 2006
  !**************************************************************

      use Driver_interface, ONLY : Driver_abortFlash

      use RuntimeParameters_interface, ONLY : RuntimeParameters_get

      use IncompNS_data, ONLY : ins_iConvU

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
      REAL, INTENT(IN):: ru1, dx, dy, dz, gravX, gravY, gravZ
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, wni, tv, visc, rho1x, rho2x, rho1y, rho2y 
      REAL, DIMENSION(:,:,:), INTENT(IN):: rho1z,rho2z 
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv, rw

!!$      REAL, DIMENSION(nx+1,ny,nz), INTENT(IN):: uni
!!$      REAL, DIMENSION(nx,ny+1,nz), INTENT(IN):: vni
!!$      REAL, DIMENSION(nx,ny,nz+1), INTENT(IN):: wni
!!$      REAL, DIMENSION(nx,ny,nz)  , INTENT(IN):: tv
!!$
!!$      REAL, DIMENSION(nx+1,ny,nz), INTENT(OUT):: ru
!!$      REAL, DIMENSION(nx,ny+1,nz), INTENT(OUT):: rv
!!$      REAL, DIMENSION(nx,ny,nz+1), INTENT(OUT):: rw

      INTEGER:: i, j, k
      REAL:: dx1, dy1, dz1, Mdens
      ! x-component variables
      REAL:: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
             uyplus, uyminus, uzplus, uzminus
      REAL:: dudxp, dudxm, dudyp, dudym, dudzp, dudzm, dvdxp, dvdxm, &
             dwdxp, dwdxm
      REAL:: tvjp, tvjm, tvkp, tvkm
      REAL:: txxp, txxm, tyyp, tyym, tzzp, tzzm
      REAL:: txyp, txym, txzp, txzm
      ! additional y-component variables
      REAL:: vyplus, vyminus, vzplus, vzminus, wyplus, wyminus
      REAL:: dvdyp, dvdym, dvdzp, dvdzm, dwdyp, dwdym
      REAL:: tvip, tvim
      REAL:: tyzp, tyzm
      ! additional z-component variables
      REAL:: wzplus, wzminus
      REAL:: dwdzp, dwdzm

      REAL:: vvip,vvim,vvjp,vvjm,vvkp,vvkm

      !KPD - 1st or 3rd order upwind
      REAL :: uu,uxp,uxm,uyp,uym,uzp,uzm,uy,uz,ududx,udvdx,udwdx
      REAL :: vv,vxp,vxm,vyp,vym,vzp,vzm,vx,vz,vdudy,vdvdy,vdwdy
      REAL :: ww,wxp,wxm,wyp,wym,wzp,wzm,wx,wy,wdudz,wdvdz,wdwdz 
      REAL :: uxpp,uxmm,uypp,uymm,uzpp,uzmm,vxpp,vxmm,vypp,vymm,vzpp,vzmm,&
              wxpp,wxmm,wypp,wymm,wzpp,wzmm
      INTEGER :: iConvU
      REAL :: rConvU

      !KPD - Choose Convection Scheme 1=1st order upwind, 2=2nd order central, 3=3rd order upwind
      !iConvU = 4
      iConvU = ins_iConvU
      !print*,"iConvU VALUE: ",iConvU

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy
      dz1 = 1.0/dz

      !++++++++++  U-COMPONENT (Variable Density)  ++++++++++
      do k = kz1,kz2
         do j = jy1,jy2
            do i = ix1,ix2+1

             !=============================================================
             !KPD - 1st or 3rd Order Upwind... ============================
             !=============================================================
             uu   = uni(i,j,k)

             uxp  = uni(i+1,j,k)
             uxm  = uni(i-1,j,k)
             uxpp = uni(i+2,j,k)
             uxmm = uni(i-2,j,k)

             uyp  = uni(i,j+1,k)
             uym  = uni(i,j-1,k)
             uypp = uni(i,j+2,k)
             uymm = uni(i,j-2,k)

             uzp  = uni(i,j,k+1)
             uzm  = uni(i,j,k-1)
             uzpp = uni(i,j,k+2)
             uzmm = uni(i,j,k-2)

             uy   = 0.25*( uni(i,j,k) + uni(i+1,j,k) + uni(i,j-1,k) + uni(i+1,j-1,k) ) 
             uz   = 0.25*( uni(i,j,k) + uni(i+1,j,k) + uni(i,j,k-1) + uni(i+1,j,k-1) )

             vv   = vni(i,j,k)

             vxp  = vni(i+1,j,k)
             vxm  = vni(i-1,j,k)
             vxpp = vni(i+2,j,k)
             vxmm = vni(i-2,j,k)

             vyp  = vni(i,j+1,k)
             vym  = vni(i,j-1,k)
             vypp = vni(i,j+2,k)
             vymm = vni(i,j-2,k)

             vzp  = vni(i,j,k+1)
             vzm  = vni(i,j,k-1)
             vzpp = vni(i,j,k+2)
             vzmm = vni(i,j,k-2)

             vx  = 0.25*( vni(i,j,k) + vni(i-1,j,k) + vni(i,j+1,k) + vni(i-1,j+1,k) )
             vz  = 0.25*( vni(i,j,k) + vni(i,j,k-1) + vni(i,j+1,k) + vni(i,j+1,k-1) )

             ww   = wni(i,j,k)

             wxp  = wni(i+1,j,k)
             wxm  = wni(i-1,j,k)
             wxpp = wni(i+2,j,k)
             wxmm = wni(i-2,j,k)

             wyp  = wni(i,j+1,k)
             wym  = wni(i,j-1,k)
             wypp = wni(i,j+2,k)
             wymm = wni(i,j-2,k)

             wzp  = wni(i,j,k+1)
             wzm  = wni(i,j,k-1)
             wzpp = wni(i,j,k+2)
             wzmm = wni(i,j,k-2)

             wx  = 0.25*( wni(i,j,k) + wni(i-1,j,k) + wni(i,j,k+1) + wni(i-1,j,k+1) )
             wy  = 0.25*( wni(i,j,k) + wni(i,j-1,k) + wni(i,j,k+1) + wni(i,j-1,k+1) )

             !=============================================================
             ! u.grad(u) = uj*dui/dxj 
             !=============================================================
             !************** X-momentum ***************************************
             if (uu .gt. 0) then
                if (iConvU .eq. 1) then
                   ududx = uu*(uu-uxm)*dx1   
                elseif (iConvU .eq. 2) then
                   ududx = uu*(3.*uu - 4.*uxm + 1.*uxmm)/2.*dx1   
                elseif (iConvU .eq. 3) then
                   ududx = uu*(2.*uxp + 3.*uu - 6.*uxm + 1.*uxmm)/6.*dx1   
                end if
             else
                if (iConvU .eq. 1) then
                   ududx = uu*(uxp-uu)*dx1   
                elseif (iConvU .eq. 2) then
                   ududx = uu*(-3.*uu + 4.*uxp - 1.*uxpp)/2.*dx1   
                elseif (iConvU .eq. 3) then
                  !ududx = uu*( 2.*uxm + 3.*uu - 6.*uxp + 1.*uxpp)/6.*dx1   
                   ududx = uu*(-2.*uxm - 3.*uu + 6.*uxp - 1.*uxpp)/6.*dx1   
                end if
             end if
             if (vx .gt. 0) then
                if (iConvU .eq. 1) then
                   vdudy = vx*(uu-uym)*dy1 
                elseif (iConvU .eq. 2) then
                   vdudy = vx*(3.*uu - 4.*uym + 1.*uymm)/2.*dy1 
                elseif (iConvU .eq. 3) then
                   vdudy = vx*(2.*uyp +3.*uu - 6.*uym + 1.*uymm)/6.*dy1 
                end if
             else
                if (iConvU .eq. 1) then
                   vdudy = vx*(uyp-uu)*dy1
                elseif (iConvU .eq. 2) then
                   vdudy = vx*(-3.*uu + 4.*uyp - 1.*uypp)/2.*dy1 
                elseif (iConvU .eq. 3) then
                  !vdudy = vx*( 2.*uym + 3.*uu - 6.*uyp + 1.*uypp)/6.*dy1 
                   vdudy = vx*(-2.*uym - 3.*uu + 6.*uyp - 1.*uypp)/6.*dy1 
                end if
             end if
             if (wx .gt. 0) then
                if (iConvU .eq. 1) then
                   wdudz = wx*(uu-uzm)*dz1 
                elseif (iConvU .eq. 2) then
                   wdudz = wx*(3.*uu - 4.*uzm + 1.*uzmm)/2.*dz1 
                elseif (iConvU .eq. 3) then
                   wdudz = wx*(2.*uzp +3.*uu - 6.*uzm + 1.*uzmm)/6.*dz1 
                end if
             else
                if (iConvU .eq. 1) then
                   wdudz = wx*(uzp-uu)*dz1 
                elseif (iConvU .eq. 2) then
                   wdudz = wx*(-3.*uu + 4.*uzp - 1.*uzpp)/2.*dz1
                elseif (iConvU .eq. 3) then
                  !wdudz = wx*( 2.*uzm + 3.*uu - 6.*uzp + 1.*uzpp)/6.*dz1
                   wdudz = wx*(-2.*uzm - 3.*uu + 6.*uzp - 1.*uzpp)/6.*dz1
                end if
             end if
             !************************************************************
             !=============================================================
             !=============================================================

               ! get velocities at 1/2 locations
               uxplus  = (uni(i+1,j  ,k  ) + uni(i  ,j  ,k  ))*0.5
               uxminus = (uni(i  ,j  ,k  ) + uni(i-1,j  ,k  ))*0.5

               uyplus  = (uni(i  ,j+1,k  ) + uni(i  ,j  ,k  ))*0.5
               uyminus = (uni(i  ,j  ,k  ) + uni(i  ,j-1,k  ))*0.5

               uzplus  = (uni(i  ,j  ,k+1) + uni(i  ,j  ,k  ))*0.5
               uzminus = (uni(i  ,j  ,k  ) + uni(i  ,j  ,k-1))*0.5

               vxplus  = (vni(i  ,j+1,k  ) + vni(i-1,j+1,k  ))*0.5
               vxminus = (vni(i  ,j  ,k  ) + vni(i-1,j  ,k  ))*0.5

               wxplus  = (wni(i  ,j  ,k+1) + wni(i-1,j  ,k+1))*0.5
               wxminus = (wni(i  ,j  ,k  ) + wni(i-1,j  ,k  ))*0.5

               ! get derivatives at 1/2 locations
               dudxp = (uni(i+1,j  ,k  ) - uni(i  ,j  ,k  ))*dx1
               dudxm = (uni(i  ,j  ,k  ) - uni(i-1,j  ,k  ))*dx1
               dudyp = (uni(i  ,j+1,k  ) - uni(i  ,j  ,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dudzp = (uni(i  ,j  ,k+1) - uni(i  ,j  ,k  ))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdxp = (vni(i  ,j+1,k  ) - vni(i-1,j+1,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dwdxp = (wni(i  ,j  ,k+1) - wni(i-1,j  ,k+1))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1

               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvjp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i-1,j+1,k  ))
               tvjm = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j  ,k+1) + tv(i-1,j  ,k+1))
               tvkm = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j,  k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))

               if (visc(i,j,kz1) .lt. 0.0 .OR. visc(i,j,kz1) .gt. 1.0) then
                  print*,"ERROR2: Bad Viscosity Value in the RHS",i,j,k,visc(i,j,k)
                  call Driver_abortFlash('ERROR: Bad RHS Viscosity')
               end if

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = visc(i,j,k)
               vvim = visc(i-1,j,k)
               vvjp = 0.25*(visc(i-1,j  ,k  ) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j+1,k  ) + visc(i-1,j+1,k  ))
               vvjm = 0.25*(visc(i-1,j-1,k  ) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvkp = 0.25*(visc(i-1,j  ,k  ) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j  ,k+1) + visc(i-1,j  ,k+1))
               vvkm = 0.25*(visc(i-1,j  ,k-1) + visc(i  ,j,  k-1) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))


               ! flux of normal total stresses, mu*dU/dXj (ru1 is 1/Re)
               txxp = (ru1*vvip + 2.0*tv(i,j,k))   *dudxp    ! KPD 2*nut, but not nu 
               txxm = (ru1*vvim + 2.0*tv(i-1,j,k)) *dudxm    ! KPD 2*nut, but not nu 
               tyyp = (ru1*vvjp + tvjp)            *dudyp
               tyym = (ru1*vvjm + tvjm)            *dudym
               tzzp = (ru1*vvkp + tvkp)            *dudzp
               tzzm = (ru1*vvkm + tvkm)            *dudzm

               ! flux of cross SGS stresses, mu*dUi/dX
               txyp = tvjp*dvdxp                             ! KPD no mol visc, turb only 
               txym = tvjm*dvdxm                             ! KPD no mol visc, turb only 
               txzp = tvkp*dwdxp                             ! KPD no mol visc, turb only 
               txzm = tvkm*dwdxm                             ! KPD no mol visc, turb only 


               !- kpd - Mixture inverse density
               Mdens = ( rho1x(i,j,k) + rho2x(i,j,k) )

               !- kpd - Choice of convection discretization...
               if (iConvU .eq. 0) then
                  rConvU = - (uxplus*uxplus - uxminus*uxminus)*dx1  &! advection term
                           - (vxplus*uyplus - vxminus*uyminus)*dy1  &
                           - (wxplus*uzplus - wxminus*uzminus)*dz1 
               elseif (iConvU .eq. 1 .OR. iConvU .eq. 3 .OR. iConvU .eq. 2) then
                  rConvU = - ududx - vdudy - wdudz
               end if


               ! calculate RHS for u-momentum
               ru(i,j,k) =                                          &                              
                          !- (uxplus*uxplus - uxminus*uxminus)*dx1  &! advection term
                          !- (vxplus*uyplus - vxminus*uyminus)*dy1  &
                          !- (wxplus*uzplus - wxminus*uzminus)*dz1  &             
                             rConvU                                 &
                           + Mdens*(txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens*(tyyp - tyym)*dy1                &
                           + Mdens*(tzzp - tzzm)*dz1                &
                           + Mdens*(txyp - txym)*dy1                &! TURBULENT cross terms
                           + Mdens*(txzp - txzm)*dz1                &
                           - gravX   

            enddo
         enddo
      enddo

      !++++++++++  V-COMPONENT  ++++++++++

      do k = kz1,kz2
         do j = jy1,jy2+1
            do i = ix1,ix2

             !=============================================================
             !KPD - 1st or 3rd Order Upwind... ============================
             !=============================================================
             uu   = uni(i,j,k)

             uxp  = uni(i+1,j,k)
             uxm  = uni(i-1,j,k)
             uxpp = uni(i+2,j,k)
             uxmm = uni(i-2,j,k)

             uyp  = uni(i,j+1,k)
             uym  = uni(i,j-1,k)
             uypp = uni(i,j+2,k)
             uymm = uni(i,j-2,k)

             uzp  = uni(i,j,k+1)
             uzm  = uni(i,j,k-1)
             uzpp = uni(i,j,k+2)
             uzmm = uni(i,j,k-2)

             uy   = 0.25*( uni(i,j,k) + uni(i+1,j,k) + uni(i,j-1,k) + uni(i+1,j-1,k) ) 
             uz   = 0.25*( uni(i,j,k) + uni(i+1,j,k) + uni(i,j,k-1) + uni(i+1,j,k-1) )

             vv   = vni(i,j,k)

             vxp  = vni(i+1,j,k)
             vxm  = vni(i-1,j,k)
             vxpp = vni(i+2,j,k)
             vxmm = vni(i-2,j,k)

             vyp  = vni(i,j+1,k)
             vym  = vni(i,j-1,k)
             vypp = vni(i,j+2,k)
             vymm = vni(i,j-2,k)

             vzp  = vni(i,j,k+1)
             vzm  = vni(i,j,k-1)
             vzpp = vni(i,j,k+2)
             vzmm = vni(i,j,k-2)

             vx  = 0.25*( vni(i,j,k) + vni(i-1,j,k) + vni(i,j+1,k) + vni(i-1,j+1,k) )
             vz  = 0.25*( vni(i,j,k) + vni(i,j,k-1) + vni(i,j+1,k) + vni(i,j+1,k-1) )

             ww   = wni(i,j,k)

             wxp  = wni(i+1,j,k)
             wxm  = wni(i-1,j,k)
             wxpp = wni(i+2,j,k)
             wxmm = wni(i-2,j,k)

             wyp  = wni(i,j+1,k)
             wym  = wni(i,j-1,k)
             wypp = wni(i,j+2,k)
             wymm = wni(i,j-2,k)

             wzp  = wni(i,j,k+1)
             wzm  = wni(i,j,k-1)
             wzpp = wni(i,j,k+2)
             wzmm = wni(i,j,k-2)

             wx  = 0.25*( wni(i,j,k) + wni(i-1,j,k) + wni(i,j,k+1) + wni(i-1,j,k+1) )
             wy  = 0.25*( wni(i,j,k) + wni(i,j-1,k) + wni(i,j,k+1) + wni(i,j-1,k+1) )

             !=============================================================
             ! u.grad(u) = uj*dui/dxj 
             !=============================================================
             !************** Y-momentum ***************************************
             if (uy .gt. 0) then
                if (iConvU .eq. 1) then
                   udvdx = uy*(vv-vxm)*dx1
                elseif (iConvU .eq. 3) then
                   udvdx = uy*(2.*vxp +3.*vv - 6.*vxm + 1.*vxmm)/6.*dx1
                end if
             else
                if (iConvU .eq. 1) then
                   udvdx = uy*(vxp-vv)*dx1
                elseif (iConvU .eq. 3) then
                  !udvdx = uy*( 2.*vxm + 3.*vv - 6.*vxp + 1.*vxpp)/6.*dx1
                   udvdx = uy*(-2.*vxm - 3.*vv + 6.*vxp - 1.*vxpp)/6.*dx1
                end if
             end if
             if (vv .gt. 0) then
                if (iConvU .eq. 1) then
                   vdvdy = vv*(vv-vym)*dy1
                elseif (iConvU .eq. 3) then
                   vdvdy = vv*(2.*vyp +3.*vv - 6.*vym + 1.*vymm)/6.*dy1
                end if
             else
                if (iConvU .eq. 1) then
                   vdvdy = vv*(vyp-vv)*dy1
                elseif (iConvU .eq. 3) then
                  !vdvdy = vv*( 2.*vym + 3.*vv - 6.*vyp + 1.*vypp)/6.*dy1
                   vdvdy = vv*(-2.*vym - 3.*vv + 6.*vyp - 1.*vypp)/6.*dy1
                end if
             end if
             if (wy .gt. 0) then
                if (iConvU .eq. 1) then
                   wdvdz = wy*(vv-vzm)*dz1
                elseif (iConvU .eq. 3) then
                   wdvdz = wy*(2.*vzp +3.*vv - 6.*vzm + 1.*vzmm)/6.*dz1
                end if
             else
                if (iConvU .eq. 1) then
                   wdvdz = wy*(vzp-vv)*dz1
                elseif (iConvU .eq. 3) then
                  !wdvdz = wy*( 2.*vzm + 3.*vv - 6.*vzp + 1.*vzpp)/6.*dz1
                   wdvdz = wy*(-2.*vzm - 3.*vv + 6.*vzp - 1.*vzpp)/6.*dz1
                end if
             end if
             !************************************************************
             !=============================================================
             !=============================================================

               ! get velocities at 1/2 locations
               vxplus  = (vni(i+1,j  ,k  ) + vni(i  ,j  ,k  ))*0.5
               vxminus = (vni(i  ,j  ,k  ) + vni(i-1,j  ,k  ))*0.5

               vyplus  = (vni(i  ,j+1,k  ) + vni(i  ,j  ,k  ))*0.5
               vyminus = (vni(i  ,j  ,k  ) + vni(i  ,j-1,k  ))*0.5

               vzplus  = (vni(i  ,j  ,k+1) + vni(i  ,j  ,k  ))*0.5
               vzminus = (vni(i  ,j  ,k  ) + vni(i  ,j  ,k-1))*0.5

               uyplus  = (uni(i+1,j  ,k  ) + uni(i+1,j-1,k  ))*0.5
               uyminus = (uni(i  ,j  ,k  ) + uni(i  ,j-1,k  ))*0.5

               wyplus  = (wni(i  ,j  ,k+1) + wni(i  ,j-1,k+1))*0.5
               wyminus = (wni(i  ,j  ,k  ) + wni(i  ,j-1,k  ))*0.5

               ! get derivatives at 1/2 locations
               dvdxp = (vni(i+1,j  ,k  ) - vni(i  ,j  ,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dvdyp = (vni(i  ,j+1,k  ) - vni(i  ,j  ,k  ))*dy1
               dvdym = (vni(i  ,j  ,k  ) - vni(i  ,j-1,k  ))*dy1
               dvdzp = (vni(i  ,j  ,k+1) - vni(i  ,j  ,k  ))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1
               dudyp = (uni(i+1,j  ,k  ) - uni(i+1,j-1,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dwdyp = (wni(i  ,j  ,k+1) - wni(i  ,j-1,k+1))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1


               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvip = 0.25*(tv(i  ,j-1,k  ) + tv(i+1,j-1,k  ) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i  ,j-1,k  ) + tv(i  ,j-1,k+1) + &
                            tv(i  ,j  ,k+1) + tv(i  ,j  ,k  ))
               tvkm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = 0.25*(visc(i  ,j-1,k  ) + visc(i+1,j-1,k  ) + &
                            visc(i+1,j  ,k  ) + visc(i  ,j  ,k  ))
               vvim = 0.25*(visc(i-1,j-1,k  ) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvjp = visc(i,j,k) 
               vvjm = visc(i,j-1,k)
               vvkp = 0.25*(visc(i  ,j-1,k  ) + visc(i  ,j-1,k+1) + &
                            visc(i  ,j  ,k+1) + visc(i  ,j  ,k  ))
               vvkm = 0.25*(visc(i  ,j-1,k-1) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i  ,j  ,k-1))

               ! flux of normal total stresses (ru1 is 1/Re)
               txxp = (ru1*vvip + tvip           )*dvdxp
               txxm = (ru1*vvim + tvim           )*dvdxm
               tyyp = (ru1*vvjp + 2.0*tv(i,j,k)  )*dvdyp
               tyym = (ru1*vvjm + 2.0*tv(i,j-1,k))*dvdym
               tzzp = (ru1*vvkp + tvkp           )*dvdzp
               tzzm = (ru1*vvkm + tvkm           )*dvdzm

               ! flux of cross SGS stresses
               txyp = tvip*dudyp
               txym = tvim*dudym
               tyzp = tvkp*dwdyp
               tyzm = tvkm*dwdym

               !- kpd - Mixture inverse density
               Mdens = ( rho1y(i,j,k) + rho2y(i,j,k) ) 

               !- kpd - Choice of convection discretization...
               if (iConvU .eq. 0) then
                  rConvU = - (uyplus*vxplus - uyminus*vxminus)*dx1   &! advection term
                           - (vyplus*vyplus - vyminus*vyminus)*dy1   &
                           - (wyplus*vzplus - wyminus*vzminus)*dz1 
               elseif (iConvU .eq. 1 .OR. iConvU .eq. 3 .OR. iConvU .eq. 2) then
                  rConvU = - udvdx - vdvdy - wdvdz
               end if

               ! calculate RHS for v-momentum
               rv(i,j,k) =                                           &
                          !- (uyplus*vxplus - uyminus*vxminus)*dx1   &! advection term
                          !- (vyplus*vyplus - vyminus*vyminus)*dy1   &
                          !- (wyplus*vzplus - wyminus*vzminus)*dz1   &
                             rConvU                                  &
                           + Mdens* (txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens* (tyyp - tyym)*dy1                &
                           + Mdens* (tzzp - tzzm)*dz1                &
                           + Mdens* (txyp - txym)*dx1                &! diffusion - cross terms
                           + Mdens* (tyzp - tyzm)*dz1                &
                           - gravY
            enddo
         enddo
      enddo

      !++++++++++  W-COMPONENT  ++++++++++
      
      do k = kz1,kz2+1
         do j = jy1,jy2
            do i = ix1,ix2

             !=============================================================
             !KPD - 1st or 3rd Order Upwind... ============================
             !=============================================================
             uu   = uni(i,j,k)

             uxp  = uni(i+1,j,k)
             uxm  = uni(i-1,j,k)
             uxpp = uni(i+2,j,k)
             uxmm = uni(i-2,j,k)

             uyp  = uni(i,j+1,k)
             uym  = uni(i,j-1,k)
             uypp = uni(i,j+2,k)
             uymm = uni(i,j-2,k)

             uzp  = uni(i,j,k+1)
             uzm  = uni(i,j,k-1)
             uzpp = uni(i,j,k+2)
             uzmm = uni(i,j,k-2)

             uy   = 0.25*( uni(i,j,k) + uni(i+1,j,k) + uni(i,j-1,k) + uni(i+1,j-1,k) ) 
             uz   = 0.25*( uni(i,j,k) + uni(i+1,j,k) + uni(i,j,k-1) + uni(i+1,j,k-1) )

             vv   = vni(i,j,k)

             vxp  = vni(i+1,j,k)
             vxm  = vni(i-1,j,k)
             vxpp = vni(i+2,j,k)
             vxmm = vni(i-2,j,k)

             vyp  = vni(i,j+1,k)
             vym  = vni(i,j-1,k)
             vypp = vni(i,j+2,k)
             vymm = vni(i,j-2,k)

             vzp  = vni(i,j,k+1)
             vzm  = vni(i,j,k-1)
             vzpp = vni(i,j,k+2)
             vzmm = vni(i,j,k-2)

             vx  = 0.25*( vni(i,j,k) + vni(i-1,j,k) + vni(i,j+1,k) + vni(i-1,j+1,k) )
             vz  = 0.25*( vni(i,j,k) + vni(i,j,k-1) + vni(i,j+1,k) + vni(i,j+1,k-1) )

             ww   = wni(i,j,k)

             wxp  = wni(i+1,j,k)
             wxm  = wni(i-1,j,k)
             wxpp = wni(i+2,j,k)
             wxmm = wni(i-2,j,k)

             wyp  = wni(i,j+1,k)
             wym  = wni(i,j-1,k)
             wypp = wni(i,j+2,k)
             wymm = wni(i,j-2,k)

             wzp  = wni(i,j,k+1)
             wzm  = wni(i,j,k-1)
             wzpp = wni(i,j,k+2)
             wzmm = wni(i,j,k-2)

             wx  = 0.25*( wni(i,j,k) + wni(i-1,j,k) + wni(i,j,k+1) + wni(i-1,j,k+1) )
             wy  = 0.25*( wni(i,j,k) + wni(i,j-1,k) + wni(i,j,k+1) + wni(i,j-1,k+1) )

             !=============================================================
             ! u.grad(u) = uj*dui/dxj 
             !=============================================================
             !************** Z-momentum ***************************************
             if (uz .gt. 0) then
                if (iConvU .eq. 1) then
                   udwdx = uz*(ww-wxm)*dz1
                elseif (iConvU .eq. 3) then
                   udwdx = uz*(2.*wxp +3.*ww - 6.*wxm + 1.*wxmm)/6.*dx1
                end if
             else
                if (iConvU .eq. 1) then
                   udwdx = uz*(wxp-ww)*dz1
                elseif (iConvU .eq. 3) then
                  !udwdx = uz*( 2.*wxm + 3.*ww - 6.*wxp + 1.*wxpp)/6.*dx1
                   udwdx = uz*(-2.*wxm - 3.*ww + 6.*wxp - 1.*wxpp)/6.*dx1
                end if
             end if
             if (vz .gt. 0) then
                if (iConvU .eq. 1) then
                   vdwdy = vz*(ww-wym)*dz1
                elseif (iConvU .eq. 3) then
                   vdwdy = vz*(2.*wyp +3.*ww - 6.*wym + 1.*wymm)/6.*dy1
                end if
             else
                if (iConvU .eq. 1) then
                   vdwdy = vz*(wyp-ww)*dz1
                elseif (iConvU .eq. 3) then
                  !vdwdy = vz*( 2.*wym + 3.*ww - 6.*wyp + 1.*wypp)/6.*dy1
                   vdwdy = vz*(-2.*wym - 3.*ww + 6.*wyp - 1.*wypp)/6.*dy1
                end if
             end if
             if (ww .gt. 0) then
                if (iConvU .eq. 1) then
                   wdwdz = ww*(ww-wzm)*dz1
                elseif (iConvU .eq. 3) then
                   wdwdz = ww*(2.*wzp +3.*ww - 6.*wzm + 1.*wzmm)/6.*dz1
                end if
             else
                if (iConvU .eq. 1) then
                   wdwdz = ww*(wzp-ww)*dz1
                elseif (iConvU .eq. 3) then
                  !wdwdz = ww*( 2.*wzm + 3.*ww - 6.*wzp + 1.*wzpp)/6.*dz1
                   wdwdz = ww*(-2.*wzm - 3.*ww + 6.*wzp - 1.*wzpp)/6.*dz1
                end if
             end if
             !************************************************************
             !=============================================================
             !=============================================================

               ! get velocities at 1/2 locations
               wxplus  = (wni(i+1,j  ,k  ) + wni(i  ,j  ,k  ))*0.5
               wxminus = (wni(i  ,j  ,k  ) + wni(i-1,j  ,k  ))*0.5
               
               wyplus  = (wni(i  ,j+1,k  ) + wni(i  ,j  ,k  ))*0.5
               wyminus = (wni(i  ,j  ,k  ) + wni(i  ,j-1,k  ))*0.5
               
               wzplus  = (wni(i  ,j  ,k+1) + wni(i  ,j  ,k  ))*0.5
               wzminus = (wni(i  ,j  ,k  ) + wni(i  ,j  ,k-1))*0.5

               uzplus  = (uni(i+1,j  ,k  ) + uni(i+1,j  ,k-1))*0.5
               uzminus = (uni(i  ,j  ,k  ) + uni(i  ,j  ,k-1))*0.5

               vzplus  = (vni(i  ,j+1,k  ) + vni(i  ,j+1,k-1))*0.5
               vzminus = (vni(i  ,j  ,k  ) + vni(i  ,j  ,k-1))*0.5

               ! get derivatives at 1/2 locations
               dwdxp = (wni(i+1,j  ,k  ) - wni(i  ,j  ,k  ))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1
               dwdyp = (wni(i  ,j+1,k  ) - wni(i  ,j  ,k  ))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1
               dwdzp = (wni(i  ,j  ,k+1) - wni(i  ,j  ,k  ))*dz1
               dwdzm = (wni(i  ,j  ,k  ) - wni(i  ,j  ,k-1))*dz1
               dudzp = (uni(i+1,j  ,k  ) - uni(i+1,j  ,k-1))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdzp = (vni(i  ,j+1,k  ) - vni(i  ,j+1,k-1))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1


               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvip = 0.25*(tv(i  ,j  ,k-1) + tv(i+1,j  ,k-1) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j  ,k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvjp = 0.25*(tv(i  ,j  ,k-1) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i  ,j+1,k-1))
               tvjm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + & 
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = 0.25*(visc(i  ,j  ,k-1) + visc(i+1,j  ,k-1) + &
                            visc(i+1,j  ,k  ) + visc(i  ,j  ,k  ))
               vvim = 0.25*(visc(i-1,j  ,k-1) + visc(i  ,j  ,k-1) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvjp = 0.25*(visc(i  ,j  ,k-1) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j+1,k  ) + visc(i  ,j+1,k-1))
               vvjm = 0.25*(visc(i  ,j-1,k-1) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i  ,j  ,k-1))
               vvkp = visc(i,j,k)
               vvkm = visc(i,j,k-1)

               ! flux of normal total stresses (ru1 is 1/Re)
               txxp = (ru1*vvip + tvip           )*dwdxp
               txxm = (ru1*vvim + tvim           )*dwdxm
               tyyp = (ru1*vvjp + tvjp           )*dwdyp
               tyym = (ru1*vvjm + tvjm           )*dwdym
               tzzp = (ru1*vvkp + 2.0*tv(i,j,k)  )*dwdzp
               tzzm = (ru1*vvkm + 2.0*tv(i,j,k-1))*dwdzm

               ! flux of cross SGS stresses
               txzp = tvip*dudzp
               txzm = tvim*dudzm
               tyzp = tvjp*dvdzp
               tyzm = tvjm*dvdzm

               !- kpd - Mixture inverse density
               Mdens = ( rho1z(i,j,k) + rho2z(i,j,k) )

               !- kpd - Choice of convection discretization...
               if (iConvU .eq. 0) then
                  rConvU = - (uzplus*wxplus - uzminus*wxminus)*dx1   &! advection term
                           - (vzplus*wyplus - vzminus*wyminus)*dy1   &
                           - (wzplus*wzplus - wzminus*wzminus)*dz1 
               elseif (iConvU .eq. 1 .OR. iConvU .eq. 3 .OR. iConvU .eq. 2) then
                  rConvU = - udwdx - vdwdy - wdwdz
               end if

               ! calculate RHS for w-momentum
               rw(i,j,k) =                                          &
                         !- (uzplus*wxplus - uzminus*wxminus)*dx1   &! advection term
                         !- (vzplus*wyplus - vzminus*wyminus)*dy1   &
                         !- (wzplus*wzplus - wzminus*wzminus)*dz1   &
                            rConvU                                  &
                          + Mdens* (txxp - txxm)*dx1                &! diffusion - normal terms
                          + Mdens* (tyyp - tyym)*dy1                &
                          + Mdens* (tzzp - tzzm)*dz1                &
                          + Mdens* (txzp - txzm)*dx1                &! diffusion - cross terms
                          + Mdens* (tyzp - tyzm)*dy1                &
                          - gravZ                   
            enddo
         enddo
      enddo

      END SUBROUTINE ins_rhs3d_VD

