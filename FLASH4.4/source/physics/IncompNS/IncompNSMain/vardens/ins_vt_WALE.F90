!!$c
!!$c-----------------------------------------------------------------------
!!$c                 ***************************                         
!!$c                 *         turvis.f        *                        
!!$c                 ***************************                       
!!$c----------------------------------------------------------------------- 
!!$c
!!$c    - Turvis:       calls the different eddy viscosity routines
!!$c
!!$c----------------------------------------------------------------------- 
!!$c
!!$c
!!$c-----SUBROUTINE-Turvis------------------------E. Balaras 7/12/98-------
!!$c
!!$c-----Adapted to AMR: M. Vanella, June 2007.----------------------------
!!$c
!!$c
  SUBROUTINE ins_vt_WALE(isgs,ng,nxc,nyc,nzc,RU1,dx,dy,dz,   &
                    coord,bsize,&
                    facexData,&
                    faceyData,&
                    facezData,&
                    solnData)             

  use ins_interface, only : ins_velgradtensor


  implicit none
 
#include "constants.h"
#include "Flash.h"


  integer isgs,ng,nxc,nyc,nzc
  real  RU1,dx,dy,dz       
  real coord(3),bsize(3)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  ! Local variables:
  real  dxdydz

  !KPD - WALE Variables...
  !=======================
  integer i,j,k,ii,jj,kk,ll
  real gkk_2,Sij_2,Sijd_2,Cm2, Oij_2, IVSO
  real, dimension(3,3) :: delta, Sij, Sijd, Oij, Otest, Stest
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: dudx,dudy,dudz, &
        dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,3,3) :: g, gij_2
  !========================================================================


  real ALPHA,RPLUS

  real ywall,ypos,dyaux
  real, parameter :: CSS = 0.16
  integer nxi,nyj,nzk

!!$c
!!$c----------------------------------------------------------------------
!!$c                      square of the ratio between test and grid filter
!!$c----------------------------------------------------------------------

! This subroutine will cause a problem for 2D flows.
! Shizhao turned off this function for 2D flows
#if NDIM == 3   
  nxi = NXB + ng
  nyj = NYB + ng
  nzk = NZB + ng


  dxdydz = (dx*dy*dz)**(0.33333333333333333)
!!$c
!!$c----------------------------------------------------------------------
!!$c                       calculate Strain Rates and |S| at cell centers
!!$c----------------------------------------------------------------------

!  call ins_velgradtensor(ng,facexData,faceyData,facezData, &              
!                         dx,dy,dz,SXX,SXY,SXZ,             &
!                         tpdvdxc,SYY,SYZ,tpdwdxc,G,SZZ)

  call ins_velgradtensor(ng,facexData,faceyData,facezData, &              
                         dx,dy,dz,dudx,dudy,dudz,             &
                         dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)

  !KPD - g(:,:,i,j) = du_i / dx_j
  !==============================
  g(:,:,:,1,1) = dudx
  g(:,:,:,1,2) = dudy
  g(:,:,:,1,3) = dudz
  g(:,:,:,2,1) = dvdx
  g(:,:,:,2,2) = dvdy
  g(:,:,:,2,3) = dvdz
  g(:,:,:,3,1) = dwdx
  g(:,:,:,3,2) = dwdy
  g(:,:,:,3,3) = dwdz

do k=1,3
  print*,"G Tensor",g(8,8,8,k,1),g(8,8,8,k,2),g(8,8,8,k,3)
end do

  delta = 0.0
  gij_2 = 0.0  
  gkk_2 = 0.0  

  delta(1,1) = 1.0
  delta(2,2) = 1.0
  delta(3,3) = 1.0

  Cm2 = 10.6*(0.15**2.0)

  do k=GRID_KLO_GC,GRID_KHI_GC
     do j=GRID_JLO_GC,GRID_JHI_GC
        do i=GRID_ILO_GC,GRID_IHI_GC

        !KPD - Initialize to 0.0 for current cell
        !========================================
        gkk_2       = 0.
        Sijd(:,:)   = 0.
        Sij (:,:)   = 0.

        !KPD - gkk_2 = g_ii*g_ii = Sum of Digonals 
        !================================================
        !gkk_2        = ( g(i,j,k,1,1) + &
        !                 g(i,j,k,2,2) + &
        !                 g(i,j,k,3,3) )**2.0
        gkk_2        = g(i,j,k,1,1)*g(i,j,k,1,1) + &
                       g(i,j,k,2,2)*g(i,j,k,2,2) + &
                       g(i,j,k,3,3)*g(i,j,k,3,3) 

        !if (k.eq.8 .AND. j.eq.8 .AND. i.eq.8) then
        !   print*,"gkk_2",gkk_2
        !end if

        !KPD - Sijd = S_ij_d = 1/2*((g_ij)^2+(g_ji)^2) - 1/3*delta_ij*(g_kk)^2
        !=====================================================================
        !!One of two ways to compute Sijd^2...
        !Sijd(:,:) = 1./2.*( (MATMUL(g(i,j,k,:,:),g(i,j,k,:,:))) + & 
        !                    (TRANSPOSE(MATMUL(g(i,j,k,:,:),g(i,j,k,:,:)))) ) &
        !          - 1./3.*( delta*gkk_2 )

        !KPD - Sij = S_ij = 1/2 * (du_i/dx_j + du_j/dx_i)
        !================================================
        Sij (:,:) = 0.5*(g(i,j,k,:,:) + TRANSPOSE(g(i,j,k,:,:)))

        !KPD - Oij = O_ij = 1/2 * (du_i/dx_j - du_j/dx_i)
        !================================================
        Oij (:,:) = 0.5*(g(i,j,k,:,:) - TRANSPOSE(g(i,j,k,:,:)))

        !if (k.eq.8 .AND. j.eq.8 .AND. i.eq.8) then
        !   do kk=1,3
        !     print*,"Sij Tensor",Sij(kk,1),Sij(kk,2),Sij(kk,3)
        !   end do
        !end if
        !if (k.eq.8 .AND. j.eq.8 .AND. i.eq.8) then
        !   do kk=1,3
        !     print*,"Oij Tensor",Oij(kk,1),Oij(kk,2),Oij(kk,3)
        !   end do
        !end if

        !KPD - Sij_2  = Sij*Sij
        !=====================
        !KPD - Sijd_2 = Sij_d*Sij_d
        !=======================
        Sij_2  = 0.
        Sijd_2 = 0.
        Oij_2  = 0.
        do ii=1,3
           do jj=1,3

              !!One of two ways to compute Sijd^2...
              !Sijd_2 = Sijd_2 + Sijd(ii,jj)*Sijd(ii,jj)
              
              Sij_2  = Sij_2  + Sij (ii,jj)*Sij (ii,jj)

              Oij_2  = Oij_2  + Oij (ii,jj)*Oij (ii,jj)

              !!Stest = Sik * Skj
              !Stest(ii,jj) = Sij(jj,1)*Sij(1,ii)+Sij(jj,2)*Sij(2,ii)+Sij(jj,3)*Sij(3,ii)
              !!Otest = Ojk * Oki
              !Otest(ii,jj) = Oij(jj,1)*Oij(1,ii)+Oij(jj,2)*Oij(2,ii)+Oij(jj,3)*Oij(3,ii)

           end do
        end do
        !if (k.eq.8 .AND. j.eq.8 .AND. i.eq.8) then
        !do kk=1,3
        !  print*,"Stest Tensor",Stest(kk,1),Stest(kk,2),Stest(kk,3)
        !end do
        !end if
        !if (k.eq.8 .AND. j.eq.8 .AND. i.eq.8) then
        !do kk=1,3
        !  print*,"Otest Tensor",Otest(kk,1),Otest(kk,2),Otest(kk,3)
        !end do
        !end if

        IVSO = 0.
        do ii=1,3
           do jj=1,3
              do kk=1,3
                 do ll=1,3

                     IVSO = IVSO + Sij(ii,kk)*Sij(kk,jj)*Oij(jj,ll)*Oij(ll,ii) 

                 end do
              end do
           end do
        end do

        Sijd_2 = 1./6.*(Sij_2*Sij_2+Oij_2*Oij_2) + & 
                 2./3.*(Sij_2*Oij_2)             + &
                 2.*IVSO

        !if (k.eq.8 .AND. j.eq.8 .AND. i.eq.8) then
        !   print*,"Sij^2=",Sij_2,"Oij^2=",Oij_2,"IVSO=",IVSO
        !   print*,"Sijd^2=",Sijd_2,"Sij^2=",Sij_2,"V=",dxdydz
        !end if

        !KPD - Turbulent Viscosity for current cell
        !==========================================
        solnData(TVIS_VAR,i,j,k) = (Cm2*dxdydz**2.0) *         & 
                                    ( Sijd_2**(3./2.) /   &
                                     (Sij_2 **(5./2.) +   & 
                                      Sijd_2**(5./4.)) ) 

        if (k.eq.8 .AND. j.eq.8 .AND. i.eq.8) then
           print*,"The Final NuT = ",solnData(TVIS_VAR,i,j,k)
        end if
 
        end do
     end do
  end do


!  ! Smagorinsky model
!  IF(ISGS==1) THEN
!
!     DO j= ng+1 , nyj
!
!        !KPD - NU_sgs = (C_s*Delta)^2 * MAG(Sij)
!        solnData(TVIS_VAR,ng+1:nxi,j,ng+1:nzk) = &
!        (CSS*DXDYDZ)**2*G(ng+1:nxi,j,ng+1:nzk)
!
!     ENDDO
!
!  endif

#endif

  RETURN

END SUBROUTINE ins_VT_WALE


