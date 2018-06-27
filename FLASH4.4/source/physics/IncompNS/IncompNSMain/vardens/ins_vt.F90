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
  SUBROUTINE ins_vt(isgs,ng,nxc,nyc,nzc,RU1,dx,dy,dz,   &
                    coord,bsize,&
                    facexData,&
                    faceyData,&
                    facezData,&
                    solnData)             

  use ins_interface, only : ins_velgradtensor

  use Multiphase_data, only: mph_rho1,mph_rho2

  implicit none
 
#include "constants.h"
#include "Flash.h"


  integer isgs,ng,nxc,nyc,nzc
  real  RU1,dx,dy,dz       
  real coord(3),bsize(3)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  ! Local variables:
  real  dxdydz
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: UC,VC,WC
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: SXX,SYY,SZZ,&
        SXY,SYZ,SXZ,tpdvdxc,tpdwdxc,G

  integer i,j,k

  real ALPHA,RPLUS

  real ywall,ypos,dyaux
  real, parameter :: css = 0.16
  integer nxi,nyj,nzk

  real :: sgn, rhoijk
!!$c
!!$c----------------------------------------------------------------------
!!$c                      square of the ratio between test and grid filter
!!$c----------------------------------------------------------------------

  nxi = NXB + ng
  nyj = NYB + ng
  nzk = NZB + ng


  dxdydz = (dx*dy*dz)**(0.33333333333333333)
!!$c
!!$c----------------------------------------------------------------------
!!$c                       calculate Strain Rates and |S| at cell centers
!!$c----------------------------------------------------------------------

!du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz, dw/dx, dw/dy, dw/dz 
  call ins_velgradtensor(ng,facexData,faceyData,facezData, &              
                         dx,dy,dz,SXX,SXY,SXZ,             &
                         tpdvdxc,SYY,SYZ,tpdwdxc,G,SZZ)

  SXY = 0.5*(SXY + tpdvdxc)   !.5*(du/dy+dv/dx) = .5*(S12+S21)
  SXZ = 0.5*(SXZ + tpdwdxc)   !.5*(du/dz+dw/dx) = .5*(S13+S31)
  SYZ = 0.5*(SYZ + G)         !.5*(dv/dz+dw/dy) = .5*(S23+S32)

  ! = SQRT( 2*SijSij )
  G = SQRT(2.*(SXX**2+SYY**2+SZZ**2)+4.*(SXY**2+SYZ**2+SXZ**2)) 


  ! Smagorinsky model
  IF(ISGS==1) THEN

     !DO j= ng+1 , nyj
     !   solnData(TVIS_VAR,ng+1:nxi,j,ng+1:nzk) = &
     !   (CSS*DXDYDZ)**2*G(ng+1:nxi,j,ng+1:nzk)
     !ENDDO

     do k=GRID_KLO,GRID_KHI
        do j=GRID_JLO,GRID_JHI
           do i=GRID_ILO,GRID_IHI

              ! Density of cell:
              sgn = sign(1.,solnData(DFUN_VAR,i,j,k))
              rhoijk = 0.5*( (1.+sgn)*mph_rho1 + (1.-sgn)*mph_rho2 )

              solnData(TVIS_VAR,i,j,k) = rhoijk*(CSS*DXDYDZ)**2*G(i,j,k)


           enddo
        enddo
     enddo



  endif



  RETURN

END SUBROUTINE ins_VT


