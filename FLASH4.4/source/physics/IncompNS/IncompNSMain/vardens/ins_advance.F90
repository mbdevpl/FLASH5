SUBROUTINE ins_predictor(uni,vni,wni,unew,vnew,wnew,uold,vold,&
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

      INTEGER :: i,j,k
      
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

!do j=jy1,jy2
!do i=ix1,ix2
!print*,"PRED",uni(i,j,1),vni(i,j,1)
!end do
!end do

#if NDIM == 3
      wni(ix1:ix2,jy1:jy2,kz1:kz2+1) = wni(ix1:ix2,jy1:jy2,kz1:kz2+1) + &
         dt*(gama*wnew(ix1:ix2,jy1:jy2,kz1:kz2+1) +                     &
             rhoa*wold(ix1:ix2,jy1:jy2,kz1:kz2+1) -                     &
             ins_prescoeff*alfa*( p(ix1:ix2,jy1:jy2,kz1:kz2+1) -        &
                                  p(ix1:ix2,jy1:jy2,kz1-1:kz2) )/dz)
#endif

END SUBROUTINE ins_predictor

!########################################################################

SUBROUTINE ins_divergence(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
         dx,dy,dz,divv)

      ! This routine computes the divergence of the velocity field.

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv

      INTEGER :: i,j,k

      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
         ( uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
           uni(ix1:ix2,jy1:jy2,kz1:kz2) )/dx +  &       
         ( vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
           vni(ix1:ix2,jy1:jy2,kz1:kz2) )/dy

!do j=jy1,jy2
!do i=ix1,ix2
!print*,"DIV",i,j,divv(i,j,1)
!end do
!end do

#if NDIM == 3
      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
           divv(ix1:ix2,jy1:jy2,kz1:kz2)     +  &
         ( wni(ix1:ix2,jy1:jy2,kz1+1:kz2+1) -   &
           wni(ix1:ix2,jy1:jy2,kz1:kz2) )/dz 
#endif

END SUBROUTINE ins_divergence

!########################################################################

SUBROUTINE ins_corrector(uni,vni,wni,p,ix1,ix2,jy1,jy2,kz1,kz2, &
        dt,dx,dy,dz,alfa)

      ! This routine computes the corrected divergence-free velocities.
    
      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz,alfa
      REAL, DIMENSION(:,:,:), INTENT(IN) :: p
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni

      REAL :: coef

      coef = dt*alfa

      uni(ix1+1:ix2,jy1:jy2,kz1:kz2) =             & 
         uni(ix1+1:ix2,jy1:jy2,kz1:kz2) -          &
         coef*( p(ix1+1:ix2,jy1:jy2,kz1:kz2) -     &
                p(ix1:ix2-1,jy1:jy2,kz1:kz2) )/dx

      vni(ix1:ix2,jy1+1:jy2,kz1:kz2) =             &
         vni(ix1:ix2,jy1+1:jy2,kz1:kz2) -          &
         coef*( p(ix1:ix2,jy1+1:jy2,kz1:kz2) -     &
                p(ix1:ix2,jy1:jy2-1,kz1:kz2) )/dy

#if NDIM == 3
      wni(ix1:ix2,jy1:jy2,kz1+1:kz2) =             &
         wni(ix1:ix2,jy1:jy2,kz1+1:kz2) -          &
         coef*( p(ix1:ix2,jy1:jy2,kz1+1:kz2) -     &
                p(ix1:ix2,jy1:jy2,kz1:kz2-1) )/dz
#endif

END SUBROUTINE ins_corrector

