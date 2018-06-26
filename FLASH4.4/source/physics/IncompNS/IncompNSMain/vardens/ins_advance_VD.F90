SUBROUTINE ins_predictor_VD(uni,vni,wni,unew,vnew,wnew,uold,vold,&
        wold,p,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,gama,rhoa,alfa)

      ! This routine computes the intermediate velocities based on
      ! the explicit second-order Adams-Bashforth scheme (gamma=1.5,
      ! rhoa=-0.5,alfa=1), or a third order Runge-Kutta method.

      use IncompNS_data, ONLY : ins_prescoeff

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: unew,vnew,wnew,&
                                            uold,vold,wold,&
                                            p
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni

      REAL :: gama,rhoa,alfa

      !- kpd - Output of pressure solution method to the screen
      !if (ins_prescoeff .gt. 0.) then
      !   print*,"KPD - Pressure Correction scheme IS in use (ins_advance.F90)."
      !else
      !   print*,"KPD - Pressure Correction scheme is NOT in use (ins_advance.F90)."
      !end if

      !---------------------------------------------------------------------
      !---------------------------------------------------------------------
      !- kpd - When pressure correction scheme is not used ins_prescoeff = 0
      !           and dP/dx is not used in u* predictor step
      !        uni and vni are u* and v*
      !---------------------------------------------------------------------
      !---------------------------------------------------------------------

      uni(ix1:ix2+1,jy1:jy2,kz1:kz2) = uni(ix1:ix2+1,jy1:jy2,kz1:kz2) + &
         dt*(gama*unew(ix1:ix2+1,jy1:jy2,kz1:kz2) +                     &
             rhoa*uold(ix1:ix2+1,jy1:jy2,kz1:kz2) -                     &
             ins_prescoeff*alfa*( p(ix1:ix2+1,jy1:jy2,kz1:kz2) -        &
                                  p(ix1-1:ix2,jy1:jy2,kz1:kz2) )/dx)

      vni(ix1:ix2,jy1:jy2+1,kz1:kz2) = vni(ix1:ix2,jy1:jy2+1,kz1:kz2) + &
         dt*(gama*vnew(ix1:ix2,jy1:jy2+1,kz1:kz2) +                     &
             rhoa*vold(ix1:ix2,jy1:jy2+1,kz1:kz2) -                     &
             ins_prescoeff*alfa*( p(ix1:ix2,jy1:jy2+1,kz1:kz2) -        &
                                  p(ix1:ix2,jy1-1:jy2,kz1:kz2) )/dy)

#if NDIM == 3
      wni(ix1:ix2,jy1:jy2,kz1:kz2+1) = wni(ix1:ix2,jy1:jy2,kz1:kz2+1) + &
         dt*(gama*wnew(ix1:ix2,jy1:jy2,kz1:kz2+1) +                     &
             rhoa*wold(ix1:ix2,jy1:jy2,kz1:kz2+1) -                     &
             ins_prescoeff*alfa*( p(ix1:ix2,jy1:jy2,kz1:kz2+1) -        &
                                  p(ix1:ix2,jy1:jy2,kz1-1:kz2) )/dz)
#endif


END SUBROUTINE ins_predictor_VD

!########################################################################

SUBROUTINE ins_divergence_VD(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
         dx,dy,dz,divv)

      ! This routine computes the divergence of the velocity field.

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv


      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
         ( uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
           uni(ix1:ix2,jy1:jy2,kz1:kz2) )/dx +  &       
         ( vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
           vni(ix1:ix2,jy1:jy2,kz1:kz2) )/dy


#if NDIM == 3
      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
           divv(ix1:ix2,jy1:jy2,kz1:kz2)     +  &
         ( wni(ix1:ix2,jy1:jy2,kz1+1:kz2+1) -   &
           wni(ix1:ix2,jy1:jy2,kz1:kz2) )/dz 
#endif

END SUBROUTINE ins_divergence_VD

!########################################################################

SUBROUTINE ins_corrector_VD(uni,vni,wni,sigx,sigy,sigz,p,ix1,ix2,jy1,jy2,kz1,kz2, &
        dt,dx,dy,dz,alfa,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z)

      ! This routine computes the corrected divergence-free velocities.
    
      use Grid_data,        ONLY : gr_meshMe

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz,alfa
      REAL, DIMENSION(:,:,:), INTENT(IN) :: p,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni,sigx,sigy,sigz

      REAL :: coef,Mdens
      INTEGER :: i,j,k

      coef = dt*alfa

      !- kpd - Doesn't loop through boundary X-mom locations.
      !        Those were taken care of in ab2rk3
      do k=kz1,kz2
         do j=jy1,jy2
            do i=ix1+1,ix2

               !- kpd - inverse of the mixture density
               !--------------------------------------
               Mdens = (rho1x(i,j,k) + rho2x(i,j,k))
               !Mdens = 1.0 !(rho1x(i,j,k) + rho2x(i,j,k))
 
               !-------------------------------------------------
               !- kpd - Correct the final velocity...
               !           u(n+1) = u(*) - dt/rho*[ dP/dx-sig*K ]
               !        
               !    *** Note: sigx already contains 1/rho(mix)
               !-------------------------------------------------
               uni(i,j,k) =                                & 
                           uni(i,j,k) -                    &
                           coef*( Mdens*( p(i,j,k) -       &
                                          p(i-1,j,k) )/dx  &
                                  -sigx(i,j,k) )

               !if (gr_meshMe .eq. 0) print*,"KPD CorrX1",i,j,k,Mdens,sigx(i,j,k),uni(i,j,k)

            enddo
         enddo
      enddo

!      uni(ix1+1:ix2,jy1:jy2,kz1:kz2) =             & 
!         uni(ix1+1:ix2,jy1:jy2,kz1:kz2) -          &
!         coef*( ( p(ix1+1:ix2,jy1:jy2,kz1:kz2) -     &
!                  p(ix1:ix2-1,jy1:jy2,kz1:kz2) )/dx  &
!               -sigx(ix1+1:ix2,jy1:jy2,kz1:kz2))


      !- kpd - Doesn't loop through boundary Y-mom locations.
      !        Those were taken care of in ab2rk3
      do k=kz1,kz2
         do j=jy1+1,jy2
            do i=ix1,ix2

               !- kpd - inverse of the mixture density
               !--------------------------------------
               Mdens = (rho1y(i,j,k) + rho2y(i,j,k))
               !Mdens = 1.0 !(rho1y(i,j,k) + rho2y(i,j,k))

               !print*,"KPD CorrY",i,j,k,Mdens,sigy(i,j,k)

               !-------------------------------------------------
               !- kpd - Correct the final velocity...
               !           u(n+1) = u(*) - dt/rho*[ dP/dx-sig*K ]
               !        
               !        Note: sigy already contains 1/rho(mix)
               !-------------------------------------------------
!print*,"CORRECT",i,j,vni(i,j,k),Mdens,(p(i,j,k)-p(i,j-1,k) )/dy,coef*Mdens*(p(i,j,k)-p(i,j-1,k) )/dy
               vni(i,j,k) =                                &           
                           vni(i,j,k) -                    &
                           coef*( Mdens*( p(i,j,k) -       &
                                          p(i,j-1,k) )/dy  &
                                  -sigy(i,j,k) )

               !if (gr_meshMe .eq. 0) print*,"KPD CorrY1",i,j,k,Mdens,sigy(i,j,k),vni(i,j,k)
            enddo
         enddo
      enddo

!      vni(ix1:ix2,jy1+1:jy2,kz1:kz2) =             &
!         vni(ix1:ix2,jy1+1:jy2,kz1:kz2) -          &
!         coef*( ( p(ix1:ix2,jy1+1:jy2,kz1:kz2) -     &
!                  p(ix1:ix2,jy1:jy2-1,kz1:kz2) )/dy  &
!               -sigy(ix1:ix2,jy1+1:jy2,kz1:kz2))

#if NDIM == 3
!      wni(ix1:ix2,jy1:jy2,kz1+1:kz2) =             &
!         wni(ix1:ix2,jy1:jy2,kz1+1:kz2) -          &
!         coef*( ( p(ix1:ix2,jy1:jy2,kz1+1:kz2) -     &
!                  p(ix1:ix2,jy1:jy2,kz1:kz2-1) )/dz  &
!               -sigz(ix1:ix2,jy1:jy2,kz1+1:kz2))

      do k=kz1+1,kz2
         do j=jy1,jy2
            do i=ix1,ix2

               !- kpd - inverse of the mixture density
               Mdens = (rho1z(i,j,k) + rho2z(i,j,k))
               !Mdens = 1.0 !(rho1z(i,j,k) + rho2z(i,j,k))

               !-------------------------------------------------
               !- kpd - Correct the final velocity...
               !           u(n+1) = u(*) - dt/rho*[ dP/dx-sig*K ]
               !        
               !        Note: sigz already contains 1/rho(mix)
               !-------------------------------------------------
               wni(i,j,k) =                                &           
                           wni(i,j,k) -                    &
                           coef*( Mdens*( p(i,j,k) -       &
                                          p(i,j,k-1) )/dz  &
                           -sigz(i,j,k))
            enddo
         enddo
      enddo

#endif

END SUBROUTINE ins_corrector_VD

