
      subroutine mph_KPDlsRedistance(s,u,v,lsit,dx,dy,ix1,ix2,jy1,jy2)

        use IncompNS_data, ONLY : ins_cfl

        implicit none

#include "Flash.h"

        !- kpd - Imported variables

        integer, intent(in) :: lsit,ix1,ix2,jy1,jy2
        real, dimension(:,:,:), intent(inout):: s
        real, dimension(:,:,:), intent(in):: u,v
        real, intent(in) :: dx,dy

        !- kpd - Local variables
        real :: so(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
                   sgn(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
                   eps,t,agf,dtL,err1,err2, &
                   s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
                   rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
                   aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
                   a1r,a2r,a3r,a1l,a2l,a3l, &
                   fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
                   frx,flx,fry,fly

        integer :: i,j,k,incrm,m,itr


        !- kpd - For 2-D simulations
        k=1

        dtL = ins_cfl / (1./dx + 1./dy)
        t = 0.0
        eps = 1E-14

        incrm = 1
        m = 1
        itr = 0

       !- kpd - Level Set Redistance Iterations (pseudo-time)
        do while(itr.lt.lsit)

           err1 = 0.

           itr = itr + 1

           print*,"Level Set Redist Iter # ",itr

           do j = jy1,jy2
              do i = ix1,ix2
                 sgn(i,j,k) = s(i,j,k)/abs(s(i,j,k)+eps)
              end do
           end do

        so = s

        !- kpd - Loop through interior cells ony within the 
        !           5th Order WENO stancil size (3 cells each way).
        !           FLASH has guard cells available.

        do j = jy1,jy2 
           do i = ix1,ix2 

              s(i,j,k) = 0.d0

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
              !- kpd - WENO3 interpolated PHI values at cell face... 
              !---------------------------------------------------------
              frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
              flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
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


              !- kpd - Compute the magnitude of the phi gradient
              !---------------------------------------------------------
              agf = sqrt( ((frx-flx)/dx)**2. + ((fry-fly)/dy)**2. )

              !---------------------------------------------------------
              !- kpd - This can be used as a level set test...
              !print*,"Magnitude of Phi Gradient Should be =1, : ",i,j,agf
              !---------------------------------------------------------

              !---------------------------------------------------------
              !- kpd - Solve Level Set distance equation ---------------
              !---------------------------------------------------------
              s(i,j,k) = so(i,j,k) + dtL*sgn(i,j,k)*(1. - agf)
              !---------------------------------------------------------

              err1 = max(err1,abs(1.-agf))

           end do
        end do

           err2 = sqrt(sum((s - so)**2)/dble(ix2-1)/dble(jy2-1))

           t = t + dtL

           if(itr.lt.lsit) then
              !write(6,'(i8,f10.5,2e20.10)') itr,t,err1
              !print*,"LS Redistance Errors: ",itr,t,err1
              !m = m + incrm
           end if

        end do


      end subroutine mph_KPDlsRedistance 



