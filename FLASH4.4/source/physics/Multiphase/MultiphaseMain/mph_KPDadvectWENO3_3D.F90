
   
     subroutine mph_KPDadvectWENO3_3D(s,u,v,w,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,blockID)

        use Simulation_data, ONLY : sim_xMax, sim_xMin, sim_yMin, &
                                    sim_yMax, sim_zMin, sim_zMax

        use RuntimeParameters_interface, ONLY : RuntimeParameters_get

        use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas
        use Grid_data, ONLY : gr_meshMe

        use IncompNS_data, ONLY : ins_gravX, ins_gravY, ins_gravZ, &
                                  ins_dampC, ins_xDampL, ins_xDampR, ins_yDampL, ins_zDampL
 
        use Driver_data, ONLY : dr_nstep

        implicit none

#include "Flash.h"
#include "constants.h"

        real, dimension(:,:,:), intent(inout):: s,u,v,w
        real, intent(in) :: dt,dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2,blockID

        real :: so(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
                   err,eps,d0,d1,d2,eyl,eyr,ur,ul,vr,vl, &
                   s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
                   rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
                   aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
                   a1r,a2r,a3r,a1l,a2l,a3l, &
                   fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
                   frx,flx,fry,fly, &
                   frz,flz,wr,wl

        integer :: i,j,k,n

        !****************************************************
        !kpd - Damping Variables...
        real :: xcell, ycell, zcell, Fn, pi, xd, xdM, xdP, ydP, ydM, zdP, zdM, &
                AA, AAx, AAy, AAz
        real del(MDIM),bsize(MDIM),coord(MDIM)
        real, dimension(2,MDIM) :: boundBox

        call RuntimeParameters_get('xmin',    sim_xMin)
        call RuntimeParameters_get('xmax',    sim_xMax)
        call RuntimeParameters_get('ymin',    sim_yMin)
        call RuntimeParameters_get('ymax',    sim_yMax)
        call RuntimeParameters_get('zmin',    sim_zMin)
        call RuntimeParameters_get('zmax',    sim_zMax)
        call RuntimeParameters_get("gravX",ins_gravX)
        call RuntimeParameters_get("gravY",ins_gravY)
        call RuntimeParameters_get("gravZ",ins_gravZ)
        call RuntimeParameters_get("dampC",ins_dampC)
        call RuntimeParameters_get("xDampL",ins_xDampL)
        call RuntimeParameters_get("xDampR",ins_xDampR)
        call RuntimeParameters_get("yDampL",ins_yDampL)
        call RuntimeParameters_get("zDampL",ins_zDampL)

        !kpd - Froude Number...
        Fn = 1./(ins_gravX+ins_gravY+ins_gravZ)

        pi = 3.14159265359

        !- kpd - Froude base damping distance...
        !xd  = sim_xMax - (2.*pi*(Fn**2.))
        !xd  = sim_xMax - 3.0 
        xd  = ins_xDampR 
        xdM  = ins_xDampL 
        xdP  = ins_xDampR 
        ydP = ins_yDampL
        ydM = -1.0*ins_yDampL
        zdP = ins_zDampL
        zdM = -1.0*ins_zDampL

        !if (gr_meshMe .eq. 0 .AND. dr_nstep .eq. 1) then
        !   print*,"Xdamp",xdM,xdP
        !   print*,"Ydamp",ydM,ydP
        !   print*,"Zdamp",zdM,zdP
        !end if

        call Grid_getDeltas(blockID,del)
        call Grid_getBlkCenterCoords(blockId,coord)
        call Grid_getBlkBoundBox(blockId,boundBox)

        bsize(:) = boundBox(2,:) - boundBox(1,:)
        !****************************************************

        eps = 1E-12

        so = s

        !- kpd - Loop through interior cells ony within the 
        !           5th Order WENO stancil size (3 cells each way).
        !           FLASH has guard cells available.

         do k = kz1,kz2 
           do j = jy1,jy2 
              do i = ix1,ix2 

              !***************** KPD *****************************************
              xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                      real(i - NGUARD - 1)*del(IAXIS) +   &
                      0.5*del(IAXIS)

              ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                      real(j - NGUARD - 1)*del(JAXIS)  +  &
                      0.5*del(JAXIS)

              zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                      real(k - NGUARD - 1)*del(KAXIS)  +  &
                      0.5*del(KAXIS)

              !***************************************************************
              !Determing damping zone...
              !***************************************************************
              AAx = 0.0
              AAy = 0.0
              AAz = 0.0

              !Outlet damping
              if (xcell .gt. xdP .AND. xdP .gt. 0.0001) then
                 AAx = ((xcell-xdP)/(sim_xMax-xdP))**2.0
              elseif (xcell .lt. xdM .AND. ABS(xdM) .gt. 0.0001) then
                 AAx = ((xcell-xdM)/(sim_xMin-xdM))**2.0
              end if

              !Spanwise Y damping
              if (ycell .gt. ydP .AND. ydP .gt. 0.0001) then
                 AAy = ((ycell-ydP)/(sim_yMax-ydP))**2.0
              elseif (ycell .lt. ydM .AND. ydM .lt. -0.0001) then
                !AAy = ((ycell-ydM)/(sim_yMax-ydM))**2.0
                 AAy = ((ycell-ydM)/(sim_yMin-ydM))**2.0
              end if

              !Spanwise Z damping
              if (zcell .gt. zdP .AND. zdP .gt. 0.0001) then
                 AAz = ((zcell-zdP)/(sim_zMax-zdP))**2.0
              elseif (zcell .lt. zdM .AND. zdM .lt. -0.0001) then
                !AAz = ((zcell-zdM)/(sim_zMax-zdM))**2.0
                 AAz = ((zcell-zdM)/(sim_zMin-zdM))**2.0
              end if

              AA = MAX(AAx,AAy,AAz)

              !***************************************************************
              !***************************************************************

              !- kpd - Velocities on faces used for divergence --> div(u*phi)
              ul = u(i,j,k)                  
              ur = u(i+1,j,k)
              vl = v(i,j,k)
              vr = v(i,j+1,k)
              wl = w(i,j,k)
              wr = w(i,j,k+1)

              !---------------------------------------------------------
              !---------------------------------------------------------
              !- kpd - WENO3 in the X-Direction ------------------------
              !---------------------------------------------------------
              !---------------------------------------------------------

              if (u(i+1,j,k) .gt. 0) then     !- kpd - u = (+) Downwind

                 s1r = so(i-2,j,k)
                 s2r = so(i-1,j,k)
                 s3r = so(i,j,k)
                 s4r = so(i+1,j,k)
                 s5r = so(i+2,j,k)   

                 rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
                       +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                 rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
                       +  1./4. *(    s2r           -    s4r )**2.
                 rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
                       +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                 aT1r = 1./10. /  ( eps + rIS1r )**2.
                 aT2r = 6./10. /  ( eps + rIS2r )**2.
                 aT3r = 3./10. /  ( eps + rIS3r )**2.

                 a1r = aT1r / ( aT1r + aT2r +aT3r )
                 a2r = aT2r / ( aT1r + aT2r +aT3r )
                 a3r = aT3r / ( aT1r + aT2r +aT3r )

                 fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
                 fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
                 fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

              else                      !- kpd - u = (-) Upwind

                 s1r = so(i-1,j,k)
                 s2r = so(i,j,k)
                 s3r = so(i+1,j,k)
                 s4r = so(i+2,j,k)
                 s5r = so(i+3,j,k)
  
                 rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
                       +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                 rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
                       +  1./4. *(    s2r           -    s4r )**2.
                 rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
                       +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                 aT1r = 3./10. /  ( eps + rIS1r )**2.
                 aT2r = 6./10. /  ( eps + rIS2r )**2.
                 aT3r = 1./10. /  ( eps + rIS3r )**2.

                 a1r = aT1r / ( aT1r + aT2r +aT3r )
                 a2r = aT2r / ( aT1r + aT2r +aT3r )
                 a3r = aT3r / ( aT1r + aT2r +aT3r )

                 fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
                 fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
                 fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

              end if

              if (u(i,j,k) .gt. 0) then     !- kpd - u = (+) Downwind  

                 s1l = so(i-3,j,k)
                 s2l = so(i-2,j,k)
                 s3l = so(i-1,j,k)
                 s4l = so(i,j,k)
                 s5l = so(i+1,j,k)   
 
                 rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
                       +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                 rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
                       +  1./4. *(    s2l           -    s4l )**2.
                 rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
                       +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                 aT1l = 1./10. /  ( eps + rIS1l )**2.
                 aT2l = 6./10. /  ( eps + rIS2l )**2.
                 aT3l = 3./10. /  ( eps + rIS3l )**2.

                 a1l = aT1l / ( aT1l + aT2l +aT3l )
                 a2l = aT2l / ( aT1l + aT2l +aT3l )
                 a3l = aT3l / ( aT1l + aT2l +aT3l )

                 fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
                 fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
                 fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

              else                      !- kpd - u = (-) Upwind

                 s1l = so(i-2,j,k)
                 s2l = so(i-1,j,k)
                 s3l = so(i,j,k)
                 s4l = so(i+1,j,k)
                 s5l = so(i+2,j,k)   

                 rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
                       +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                 rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
                       +  1./4. *(    s2l           -    s4l )**2.
                 rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
                       +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                 aT1l = 3./10. /  ( eps + rIS1l )**2.
                 aT2l = 6./10. /  ( eps + rIS2l )**2.
                 aT3l = 1./10. /  ( eps + rIS3l )**2.

                 a1l = aT1l / ( aT1l + aT2l +aT3l )
                 a2l = aT2l / ( aT1l + aT2l +aT3l )
                 a3l = aT3l / ( aT1l + aT2l +aT3l )

                 fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
                 fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
                 fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

              end if

              !---------------------------------------------------------
              !- kpd - WENO3 interpolated PHI values at cell X faces... 
              !---------------------------------------------------------
              frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
              flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
              !---------------------------------------------------------
              !---------------------------------------------------------
           


              !---------------------------------------------------------
              !---------------------------------------------------------
              !- kpd - WENO3 in the Y-Direction ------------------------
              !---------------------------------------------------------
              !---------------------------------------------------------

              if (v(i,j+1,k) .gt. 0) then     !- kpd - v = (+) Downwind

                 s1r = so(i,j-2,k)
                 s2r = so(i,j-1,k)
                 s3r = so(i,j,k)
                 s4r = so(i,j+1,k)
                 s5r = so(i,j+2,k)

                 rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
                       +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                 rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
                       +  1./4. *(    s2r           -    s4r )**2.
                 rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
                       +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                 aT1r = 1./10. /  ( eps + rIS1r )**2.
                 aT2r = 6./10. /  ( eps + rIS2r )**2.
                 aT3r = 3./10. /  ( eps + rIS3r )**2.

                 a1r = aT1r / ( aT1r + aT2r +aT3r )
                 a2r = aT2r / ( aT1r + aT2r +aT3r )
                 a3r = aT3r / ( aT1r + aT2r +aT3r )

                 fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
                 fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
                 fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

              else                      !- kpd - v = (-) Upwind

                 s1r = so(i,j-1,k)
                 s2r = so(i,j,k)
                 s3r = so(i,j+1,k)
                 s4r = so(i,j+2,k)
                 s5r = so(i,j+3,k)

                 rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
                       +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                 rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
                       +  1./4. *(    s2r           -    s4r )**2.
                 rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
                       +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                 aT1r = 3./10. /  ( eps + rIS1r )**2.
                 aT2r = 6./10. /  ( eps + rIS2r )**2.
                 aT3r = 1./10. /  ( eps + rIS3r )**2.

                 a1r = aT1r / ( aT1r + aT2r +aT3r )
                 a2r = aT2r / ( aT1r + aT2r +aT3r )
                 a3r = aT3r / ( aT1r + aT2r +aT3r )

                 fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
                 fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
                 fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

              end if

              if (v(i,j,k) .gt. 0) then     !- kpd - v = (+) Downwind

                 s1l = so(i,j-3,k)
                 s2l = so(i,j-2,k)
                 s3l = so(i,j-1,k)
                 s4l = so(i,j,k)
                 s5l = so(i,j+1,k)

                 rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
                       +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                 rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
                       +  1./4. *(    s2l           -    s4l )**2.
                 rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
                       +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                 aT1l = 1./10. /  ( eps + rIS1l )**2.
                 aT2l = 6./10. /  ( eps + rIS2l )**2.
                 aT3l = 3./10. /  ( eps + rIS3l )**2.

                 a1l = aT1l / ( aT1l + aT2l +aT3l )
                 a2l = aT2l / ( aT1l + aT2l +aT3l )
                 a3l = aT3l / ( aT1l + aT2l +aT3l )

                 fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
                 fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
                 fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

              else                      !- kpd - v = (-) Upwind

                 s1l = so(i,j-2,k)
                 s2l = so(i,j-1,k)
                 s3l = so(i,j,k)
                 s4l = so(i,j+1,k)
                 s5l = so(i,j+2,k)

                 rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
                       +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                 rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
                       +  1./4. *(    s2l           -    s4l )**2.
                 rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
                       +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                 aT1l = 3./10. /  ( eps + rIS1l )**2.
                 aT2l = 6./10. /  ( eps + rIS2l )**2.
                 aT3l = 1./10. /  ( eps + rIS3l )**2.

                 a1l = aT1l / ( aT1l + aT2l +aT3l )
                 a2l = aT2l / ( aT1l + aT2l +aT3l )
                 a3l = aT3l / ( aT1l + aT2l +aT3l )

                 fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
                 fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
                 fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

              end if

              !---------------------------------------------------------
              !- kpd - WENO3 interpolated PHI values at cell Y faces... 
              !---------------------------------------------------------
              fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
              fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
              !---------------------------------------------------------
              !---------------------------------------------------------


              !---------------------------------------------------------
              !---------------------------------------------------------
              !- kpd - WENO3 in the Z-Direction ------------------------
              !---------------------------------------------------------
              !---------------------------------------------------------

              if (w(i,j,k+1) .gt. 0) then     !- kpd - w = (+) Downwind

                 s1r = so(i,j,k-2)
                 s2r = so(i,j,k-1)
                 s3r = so(i,j,k)
                 s4r = so(i,j,k+1)
                 s5r = so(i,j,k+2)

                 rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
                       +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                 rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
                       +  1./4. *(    s2r           -    s4r )**2.
                 rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
                       +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                 aT1r = 1./10. /  ( eps + rIS1r )**2.
                 aT2r = 6./10. /  ( eps + rIS2r )**2.
                 aT3r = 3./10. /  ( eps + rIS3r )**2.

                 a1r = aT1r / ( aT1r + aT2r +aT3r )
                 a2r = aT2r / ( aT1r + aT2r +aT3r )
                 a3r = aT3r / ( aT1r + aT2r +aT3r )

                 fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
                 fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
                 fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

              else                      !- kpd - w = (-) Upwind

                 s1r = so(i,j,k-1)
                 s2r = so(i,j,k)
                 s3r = so(i,j,k+1)
                 s4r = so(i,j,k+2)
                 s5r = so(i,j,k+3)

                 rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
                       +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
                 rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
                       +  1./4. *(    s2r           -    s4r )**2.
                 rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
                       +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

                 aT1r = 3./10. /  ( eps + rIS1r )**2.
                 aT2r = 6./10. /  ( eps + rIS2r )**2.
                 aT3r = 1./10. /  ( eps + rIS3r )**2.

                 a1r = aT1r / ( aT1r + aT2r +aT3r )
                 a2r = aT2r / ( aT1r + aT2r +aT3r )
                 a3r = aT3r / ( aT1r + aT2r +aT3r )

                 fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
                 fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
                 fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

              end if

              if (w(i,j,k) .gt. 0) then     !- kpd - w = (+) Downwind

                 s1l = so(i,j,k-3)
                 s2l = so(i,j,k-2)
                 s3l = so(i,j,k-1)
                 s4l = so(i,j,k)
                 s5l = so(i,j,k+1)

                 rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
                       +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                 rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
                       +  1./4. *(    s2l           -    s4l )**2.
                 rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
                       +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                 aT1l = 1./10. /  ( eps + rIS1l )**2.
                 aT2l = 6./10. /  ( eps + rIS2l )**2.
                 aT3l = 3./10. /  ( eps + rIS3l )**2.

                 a1l = aT1l / ( aT1l + aT2l +aT3l )
                 a2l = aT2l / ( aT1l + aT2l +aT3l )
                 a3l = aT3l / ( aT1l + aT2l +aT3l )

                 fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
                 fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
                 fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

              else                      !- kpd - w = (-) Upwind

                 s1l = so(i,j,k-2)
                 s2l = so(i,j,k-1)
                 s3l = so(i,j,k)
                 s4l = so(i,j,k+1)
                 s5l = so(i,j,k+2)

                 rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
                       +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
                 rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
                       +  1./4. *(    s2l           -    s4l )**2.
                 rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
                       +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

                 aT1l = 3./10. /  ( eps + rIS1l )**2.
                 aT2l = 6./10. /  ( eps + rIS2l )**2.
                 aT3l = 1./10. /  ( eps + rIS3l )**2.

                 a1l = aT1l / ( aT1l + aT2l +aT3l )
                 a2l = aT2l / ( aT1l + aT2l +aT3l )
                 a3l = aT3l / ( aT1l + aT2l +aT3l )

                 fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
                 fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
                 fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

              end if

              !---------------------------------------------------------
              !- kpd - WENO3 interpolated PHI values at cell Z faces... 
              !---------------------------------------------------------
              frz = a1r*fT1r + a2r*fT2r + a3r*fT3r
              flz = a1l*fT1l + a2l*fT2l + a3l*fT3l
              !---------------------------------------------------------
              !---------------------------------------------------------


              !---------------------------------------------------------------
              !- kpd - Calculate the new Level Set function ------------------
              !!---------------------------------------------------------------
              !if (xcell .gt. xdP .AND. ydP .gt. 0.0001) then
              !   s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
              !                        - dt*(fry*vr - fly*vl)/dy &
              !                        - dt*(frz*wr - flz*wl)/dz &
              !                        - ins_dampC*AA*(s(i,j,k)-zcell)
              !else if (xcell .lt. xdM .AND. ydP .gt. 0.0001) then
              !   s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
              !                        - dt*(fry*vr - fly*vl)/dy &
              !                        - dt*(frz*wr - flz*wl)/dz &
              !                        - ins_dampC*AA*(s(i,j,k)-zcell)
              !else if (xcell .gt. xdP .AND. zdP .gt. 0.0001) then
              !   s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
              !                        - dt*(fry*vr - fly*vl)/dy &
              !                        - dt*(frz*wr - flz*wl)/dz &
              !                        - ins_dampC*AA*(s(i,j,k)-ycell)
              !else if (xcell .lt. xdM .AND. zdP .gt. 0.0001) then
              !   s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
              !                        - dt*(fry*vr - fly*vl)/dy &
              !                        - dt*(frz*wr - flz*wl)/dz &
              !                        - ins_dampC*AA*(s(i,j,k)-ycell)
              !else if (ABS(ycell) .gt. ydP .AND. ydP .gt. 0.0001) then
              !   s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
              !                        - dt*(fry*vr - fly*vl)/dy &
              !                        - dt*(frz*wr - flz*wl)/dz &
              !                        - ins_dampC*AA*(s(i,j,k)-zcell)
              !else if (ABS(zcell) .gt. zdP .AND. zdP .gt. 0.0001) then
              !   s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
              !                        - dt*(fry*vr - fly*vl)/dy &
              !                        - dt*(frz*wr - flz*wl)/dz &
              !                        - ins_dampC*AA*(s(i,j,k)-ycell)
              !else
              !   s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
              !                        - dt*(fry*vr - fly*vl)/dy &
              !                        - dt*(frz*wr - flz*wl)/dz 
              !end if
              !!*********************************************************
              if (xcell .gt. xdP .AND. ins_gravZ .gt. 0.0001) then
                 s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
                                      - dt*(fry*vr - fly*vl)/dy &
                                      - dt*(frz*wr - flz*wl)/dz &
                                      - ins_dampC*AA*(s(i,j,k)-zcell)
              else if (xcell .gt. xdP .AND. ins_gravY .gt. 0.0001) then
                 s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
                                      - dt*(fry*vr - fly*vl)/dy &
                                      - dt*(frz*wr - flz*wl)/dz &
                                      - ins_dampC*AA*(s(i,j,k)-ycell)
              else if (xcell .lt. xdM .AND. ins_gravZ .gt. 0.0001) then
                 s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
                                      - dt*(fry*vr - fly*vl)/dy &
                                      - dt*(frz*wr - flz*wl)/dz &
                                      - ins_dampC*AA*(s(i,j,k)-zcell)
              else if (xcell .lt. xdM .AND. ins_gravY .gt. 0.0001) then
                 s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
                                      - dt*(fry*vr - fly*vl)/dy &
                                      - dt*(frz*wr - flz*wl)/dz &
                                      - ins_dampC*AA*(s(i,j,k)-ycell)
              else
                 s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
                                      - dt*(fry*vr - fly*vl)/dy &
                                      - dt*(frz*wr - flz*wl)/dz 
              end if
              !---------------------------------------------------------------
              !---------------------------------------------------------------

              end do
           end do
        end do

        err = sqrt(sum((s - so)**2)/dble(NXB-1)/dble(NYB-1)/dble(NZB-1))


      end subroutine mph_KPDadvectWENO3_3D
