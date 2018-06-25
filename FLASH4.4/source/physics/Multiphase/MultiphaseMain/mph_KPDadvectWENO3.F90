
   
     subroutine mph_KPDadvectWENO3(s,u,v,dt,dx,dy,ix1,ix2,jy1,jy2,blockID)

        use Simulation_data, ONLY : sim_xMax, sim_yMax

        use RuntimeParameters_interface, ONLY : RuntimeParameters_get

        use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas

        use IncompNS_data, ONLY : ins_gravX, ins_gravY, ins_gravZ, &
                                  ins_dampC, ins_xDampL, ins_zDampL 

        implicit none

#include "Flash.h"
#include "constants.h"

        real, dimension(:,:,:), intent(inout):: s,u,v
        real, intent(in) :: dt,dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2,blockID

        !real :: so(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
        real :: so(NXB+2*NGUARD,NYB+2*NGUARD,1), &
                   err,eps,d0,d1,d2,eyl,eyr,ur,ul,vr,vl, &
                   s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
                   rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
                   aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
                   a1r,a2r,a3r,a1l,a2l,a3l, &
                   fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
                   frx,flx,fry,fly

        integer :: i,j,k,n

        !kpd - Damping Variables...
        real :: xcell, ycell, Fn, pi, xd, AA 
        real del(MDIM),bsize(MDIM),coord(MDIM)
        real, dimension(2,MDIM) :: boundBox

        call RuntimeParameters_get('xmax',    sim_xMax)
        call RuntimeParameters_get('ymax',    sim_yMax)
        call RuntimeParameters_get("gravX",ins_gravX)
        call RuntimeParameters_get("gravY",ins_gravY)
        call RuntimeParameters_get("gravZ",ins_gravZ)
        call RuntimeParameters_get("dampC",ins_dampC)
        call RuntimeParameters_get("xDampL",ins_xDampL)
        call RuntimeParameters_get("zDampL",ins_zDampL)

        !kpd - Froude Number...
        Fn = 1./(ins_gravX+ins_gravY+ins_gravZ)

        pi = 3.14159265359

        !- kpd - Froude base damping distance...
        !xd = sim_xMax - (2.*pi*(Fn**2.))
        xd  = ins_xDampL


        call Grid_getDeltas(blockID,del)
        call Grid_getBlkCenterCoords(blockId,coord)
        call Grid_getBlkBoundBox(blockId,boundBox)

        bsize(:) = boundBox(2,:) - boundBox(1,:)

        !- kpd - for 2D calcs
        k=1

        eps = 1E-15

        so = s

        !----------------------------------------------------------
        !- kpd - Loop through interior cells ony within the 
        !           5th Order WENO stancil size (3 cells each way).
        !           FLASH has guard cells available.
        !----------------------------------------------------------
        do j = jy1,jy2 
           do i = ix1,ix2 

              AA = 0.0

              !***************** KPD **********************
              xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                      real(i - NGUARD - 1)*del(IAXIS) +   &
                      0.5*del(IAXIS)

              ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                      real(j - NGUARD - 1)*del(JAXIS)  +  &
                      0.5*del(JAXIS)

              if (xcell .gt. xd) then
                 !- kpd - LS non-reflective damping term for boundaries 
                 AA = ((xcell-xd)/(sim_xMax-xd))**2.0
              end if
  

              !***************** KPD **********************

              !- kpd - Velocities on faces used for divergence --> div(u*phi)
              ul = u(i,j,k)                  
              ur = u(i+1,j,k)
              vl = v(i,j,k)
              vr = v(i,j+1,k)

              !---------------------------------------------------------
              !---------------------------------------------------------
              !- kpd - WENO3 in the Y-Direction ------------------------
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
              !- kpd - WENO3 interpolated PHI values at cell face... 
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

              if (v(i,j+1,k) .gt. 0) then     !- kpd - u = (+) Downwind

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

              else                      !- kpd - u = (-) Upwind

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

              if (v(i,j,k) .gt. 0) then     !- kpd - u = (+) Downwind

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

              else                      !- kpd - u = (-) Upwind

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
              !- kpd - WENO3 interpolated PHI values at cell face... 
              !---------------------------------------------------------
              fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
              fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
              !---------------------------------------------------------
              !---------------------------------------------------------


              !---------------------------------------------------------------
              !- kpd - Calculate the new Level Set function ------------------
              !---------------------------------------------------------------
              !if (blockID .eq. 15 .AND. i.eq.8 .AND. j.eq.8) then
              !   print*,"Block15 8,8:",xcell,ycell,s(i,j,k),AA*(s(i,j,k)-ycell)
              !end if
              s(i,j,k) = so(i,j,k) - dt*(frx*ur - flx*ul)/dx &
                                   - dt*(fry*vr - fly*vl)/dy &
                                   - ins_dampC*AA*(s(i,j,k)-ycell)
              !---------------------------------------------------------------
              !---------------------------------------------------------------

           end do
        end do

        err = sqrt(sum((s - so)**2)/dble(NXB-1)/dble(NYB-1))


      end subroutine mph_KPDadvectWENO3
