!=========================================================================
!=========================================================================
!=========================================================================


        subroutine mph_KPDcurvature3DAB(s,crv,dx,dy,dz, &
           ix1,ix2,jy1,jy2,kz1,kz2, &
           rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,pf,rho1,rho2,visc,vis1,vis2)

        implicit none

#include "Flash.h"

        !*****************************************************************

        !---------------------------------
        !- kpd - Data from routine call...
        !---------------------------------
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dx, dy, dz, rho1, rho2, vis1, vis2
        real, dimension(:,:,:), intent(inout):: s,crv, &
                                                rho1x,rho2x,rho1y, &
                                                rho2y,pf, &
                                                rho1z,rho2z, visc

        !--------------------------
        !- kpd - Local variables...
        !--------------------------
        integer :: i,j,k
        real :: rPhiXN,rPhiXE,rPhiXS,rPhiXW, &
                rPhiYN,rPhiYE,rPhiYS,rPhiYW, &
                rMagN,rMagE,rMagS,rMagW, &
                rPhiZN,rPhiZE,rPhiZS,rPhiZW, &
                rMagF,rMagB, &
                rPhiXF,rPhiXB,rPhiYF,rPhiYB,rPhiZF,rPhiZB

        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr, &
                zijl,zijr,zid,zij,zidl,zidr
        real, parameter :: eps = 1E-13

        !*****************************************************************

        !---------------------------------------------
        !- kpd - Compute the curvature ---------------
        !---------------------------------------------

        crv = 0.
        do k = kz1,kz2
           do j = jy1,jy2
              do i = ix1,ix2    

              !-------------------------------------------------
              !- kpd - 3 phi gradients per face method
              !-------------------------------------------------

              !- kpd - Compute [d(phi)/dx] on all faces
              rPhiXE = 1./dx*(s(i+1,j,k) - s(i,j,k)  )
              rPhiXW = 1./dx*(s(i,j,k)   - s(i-1,j,k))
              rPhiXN = 1./4./dx * ( (s(i+1,j+1,k) - s(i-1,j+1,k)) &
                                  + (s(i+1,j,k)   - s(i-1,j,k)  ) )
              rPhiXS = 1./4./dx * ( (s(i+1,j,k)   - s(i-1,j,k)  ) &
                                  + (s(i+1,j-1,k) - s(i-1,j-1,k)) )
              rPhiXF = 1./4./dx * ( (s(i+1,j,k+1) - s(i-1,j,k+1)) &
                                  + (s(i+1,j,k)   - s(i-1,j,k)  ) )
              rPhiXB = 1./4./dx * ( (s(i+1,j,k)   - s(i-1,j,k)  ) &
                                  + (s(i+1,j,k-1) - s(i-1,j,k-1)) )

              !- kpd - Compute [d(phi)/dy] on all faces
              rPhiYE = 1./4./dy * ( (s(i+1,j+1,k) - s(i+1,j-1,k)) &
                                  + (s(i,j+1,k)   - s(i,j-1,k)  ) )
              rPhiYW = 1./4./dy * ( (s(i,j+1,k)   - s(i,j-1,k)  ) &
                                  + (s(i-1,j+1,k) - s(i-1,j-1,k)) )
              rPhiYN = 1./dy*(s(i,j+1,k) - s(i,j,k)  )
              rPhiYS = 1./dy*(s(i,j,k)   - s(i,j-1,k))
              rPhiYF = 1./4./dy * ( (s(i,j+1,k+1) - s(i,j-1,k+1)) &
                                  + (s(i,j+1,k)   - s(i,j-1,k)  ) )
              rPhiYB = 1./4./dy * ( (s(i,j+1,k)   - s(i,j-1,k)  ) &
                                  + (s(i,j+1,k-1) - s(i,j-1,k-1)) )

              !- kpd - Compute [d(phi)/dz] on all faces
              rPhiZE = 1./4./dz * ( (s(i+1,j,k+1) - s(i+1,j,k-1)) &
                                  + (s(i,j,k+1)   - s(i,j,k-1)  ) )
              rPhiZW = 1./4./dz * ( (s(i,j,k+1)   - s(i,j,k-1)  ) &
                                  + (s(i-1,j,k+1) - s(i-1,j,k-1)) )
              rPhiZN = 1./4./dz * ( (s(i,j+1,k+1) - s(i,j+1,k-1)) &
                                  + (s(i,j,k+1)   - s(i,j,k-1)  ) )
              rPhiZS = 1./4./dz * ( (s(i,j,k+1)   - s(i,j,k-1)  ) &
                                  + (s(i,j-1,k+1) - s(i,j-1,k-1)) )
              rPhiZF = 1./dz*(s(i,j,k+1) - s(i,j,k)  )
              rPhiZB = 1./dz*(s(i,j,k)   - s(i,j,k-1))


              !- kpd - Compute the magnitude of the normal for ALL faces
              rMagE = sqrt( rPhiXE**2. + rPhiYE**2. + rPhiZE**2.) + eps
              rMagW = sqrt( rPhiXW**2. + rPhiYW**2. + rPhiZW**2.) + eps
              rMagN = sqrt( rPhiXN**2. + rPhiYN**2. + rPhiZN**2.) + eps
              rMagS = sqrt( rPhiXS**2. + rPhiYS**2. + rPhiZS**2.) + eps
              rMagF = sqrt( rPhiXF**2. + rPhiYF**2. + rPhiZF**2.) + eps
              rMagB = sqrt( rPhiXB**2. + rPhiYB**2. + rPhiZB**2.) + eps

              !------------------------------------------------------------
              !- kpd - Finally, compue the curvature, K=grad(s)/||grad(s)||
              crv(i,j,k) = 1./dx * (rPhiXE/rMagE - rPhiXW/rMagW) &
                         + 1./dy * (rPhiYN/rMagN - rPhiYS/rMagS) &
                         + 1./dz * (rPhiZF/rMagF - rPhiZB/rMagB) 
              !------------------------------------------------------------

              end do
           end do
        end do

        !********************************************************************

        !----------------------------------------------------
        !- kpd - Set phase function on each side of interface
        !----------------------------------------------------
        do k = kz1-1,kz2+1
           do j = jy1-1,jy2+1
              do i = ix1-1,ix2+1
                 pf(i,j,k) = 0.

                 if(s(i,j,k).ge.0.) then
                    pf(i,j,k) = 1.                       
                    visc(i,j,k) = vis1/vis2               !- kpd - Set viscosity on each side of interface
                 else
                    visc(i,j,k) = vis2/vis2
                 end if

              end do
           end do
        end do

        !********************************************************************

        !--------------------------------------------------------------
        !- kpd - These are FACE VALUED inverse densities for each phase
        !--------------------------------------------------------------

        !- kpd - density on x-faces
        !rho1x = 0.
        !rho2x = 0.
        !- kpd - Loop through boundary and interior cell faces
        do k = kz1-1,kz2+1
           do j = jy1-1,jy2+1
              do i = ix1-1,ix2+1

              rho1x(i,j,k) = 0.
              rho2x(i,j,k) = 0.

              a1 = (pf(i-1,j,k) + pf(i,j,k)) / 2.                       
              a2 = pf(i-1,j,k)  /abs(pf(i-1,j,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)
              !rho1x(i,j,k) = a1*a2/rho1
              !rho2x(i,j,k) = (1. - a1*a2)/rho2
              rho1x(i,j,k) = a1*a2/(rho1/rho2)
              rho2x(i,j,k) = (1. - a1*a2)/(rho2/rho2)

              end do
           end do
        end do

        !- kpd - density on y-faces
        !rho1y = 0.
        !rho2y = 0.
        !- kpd - Loop through boundary and interior cell faces
        do i = ix1-1,ix2+1
           do k = kz1-1,kz2+1
              do j = jy1-1,jy2+1

              rho1y(i,j,k) = 0.
              rho2y(i,j,k) = 0.

              a1 = (pf(i,j-1,k) + pf(i,j,k)) / 2.           
              a2 = pf(i,j-1,k)  /abs(pf(i,j-1,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)
              !rho1y(i,j,k) = a1*a2/rho1
              !rho2y(i,j,k) = (1. - a1*a2)/rho2
              rho1y(i,j,k) = a1*a2/(rho1/rho2)
              rho2y(i,j,k) = (1. - a1*a2)/(rho2/rho2)

              end do
           end do
        end do

        !- kpd - density on z-faces
        !rho1z = 0.
        !rho2z = 0.
        !- kpd - Loop through boundary and interior cell faces
        do i = ix1-1,ix2+1
           do j = jy1-1,jy2+1
              do k = kz1-1,kz2+1

              rho1z(i,j,k) = 0.
              rho2z(i,j,k) = 0.

              a1 = (pf(i,j,k-1) + pf(i,j,k)) / 2.           
              a2 = pf(i,j,k-1)  /abs(pf(i,j,k-1)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)
              !rho1z(i,j,k) = a1*a2/rho1
              !rho2z(i,j,k) = (1. - a1*a2)/rho2
              rho1z(i,j,k) = a1*a2/(rho1/rho2)
              rho2z(i,j,k) = (1. - a1*a2)/(rho2/rho2)

              end do
           end do
        end do

      end subroutine mph_KPDcurvature3DAB

!=========================================================================
!=========================================================================
!=========================================================================


!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!- kpd - There is a break in the ...curvature routine to 
!        go back and do guard cell filling.
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------

!=========================================================================
!=========================================================================
!=========================================================================


        subroutine mph_KPDcurvature3DC(s,crv,rho1x,rho2x,rho1y,rho2y, &
                                       pf,w,sigx,sigy,dx,dy,          &
                                       rho1,rho2,xit,ix1,ix2, &
                                       jy1,jy2,dz,kz1,kz2,rho1z, &
                                       rho2z,sigz)


        implicit none

#include "Flash.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dx, dy, dz, rho1, rho2, xit

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,sigx,sigy, &
                                                rho1z,rho2z,sigz

        !integer :: icrv(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
        integer :: icrv(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD)

        !- kpd - Local Data 
        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr, &
                zijl,zijr,zid,zij,zidl,zidr
        integer :: i,j,k
        real, parameter :: eps = 1.0E-13


        sigx = 0.
        sigy = 0.
        sigz = 0.
        w = 0.
        icrv = 0

        do k = kz1-1,kz2
           do j = jy1-1,jy2
              do i = ix1-1,ix2

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i+1,j,k).eq.1.) then

                 th = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i+1,j,k)*(1.-th) + crv(i,j,k)*th

                 xijl = xit*crv(i,j,k)                    !- kpd - sigma*K. Used for jump in pressure
                 xijr = xit*crv(i+1,j,k)                  !- kpd - sigma*K. Used for jump in pressure
                 xidl = 0.                                !- kpd - Used for jump in gradient
                 xidr = 0.                                !- kpd - Used for jump in gradient

                 xij = xijl*th + xijr*(1.-th)             !- kpd - Jump in value
                 xid = xidl*th + xidr*(1.-th)             !- kpd - Jump in gradient. Equal to 0 here.

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1x=0 and rho2x=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Mixture density (Smeared)
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 !--------------------------------------------------------------------------

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   - xij/aa/dx**2 - xid*th*(rho1/rho2)/aa/dx
                 w(i+1,j,k) = w(i+1,j,k) + xij/aa/dx**2 - xid*(1.-th)*(rho2/rho2)/aa/dx


                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = - xij/aa/dx


                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i+1,j,k).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k)*(1.-th) + crv(i+1,j,k)*th

                 xijl = xit*crv(i,j,k)
                 xijr = xit*crv(i+1,j,k)
                 xidl = 0.
                 xidr = 0.

                 xij = xijl*(1.-th) + xijr*th
                 xid = xidl*(1.-th) + xidr*th

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1x=0 and rho2x=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Density IS SMEARED HERE
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 !--------------------------------------------------------------------------


                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + xij/aa/dx**2 + xid*(1.-th)*(rho2/rho2)/aa/dx
                 w(i+1,j,k) = w(i+1,j,k) - xij/aa/dx**2 + xid*th*(rho1/rho2)/aa/dx


                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = xij/aa/dx


                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i,j+1,k).eq.1.) then

                 th = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j+1,k)*(1.-th) + crv(i,j,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0.
                 yidr = 0.

                 yij = yijl*th + yijr*(1.-th)
                 yid = yidl*th + yidr*(1.-th)

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1y=0 and rho2y=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)               !- kpd - Density IS SMEARED HERE
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa   !- kpd - Density IS SMEARED HERE
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa   !- kpd - Density IS SMEARED HERE


                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k) - yij/aa/dy**2 - yid*th*(rho1/rho2)/aa/dy
                 w(i,j+1,k) = w(i,j+1,k)   + yij/aa/dy**2 - yid*(1.-th)*(rho2/rho2)/aa/dy


                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = - yij/aa/dy


                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i,j+1,k).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k)*(1.-th) + crv(i,j+1,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0.
                 yidr = 0.

                 yij = yijl*(1.-th) + yijr*th
                 yid = yidl*(1.-th) + yidr*th

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1y=0 and rho2y=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Density IS SMEARED HERE
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE


                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + yij/aa/dy**2 + yid*(1.-th)*(rho2/rho2)/aa/dy
                 w(i,j+1,k) = w(i,j+1,k) - yij/aa/dy**2 + yid*th*(rho1/rho2)/aa/dy


                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = yij/aa/dy


                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell to front
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i,j,k+1).eq.1.) then

                 th = abs(s(i,j,k+1))/(abs(s(i,j,k+1))+abs(s(i,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k+1)*(1.-th) + crv(i,j,k)*th

                 zijl = xit*crv(i,j,k)                    !- kpd - sigma*K. Used for jump in pressure
                 zijr = xit*crv(i,j,k+1)                  !- kpd - sigma*K. Used for jump in pressure
                 zidl = 0.                                !- kpd - Used for jump in gradient
                 zidr = 0.                                !- kpd - Used for jump in gradient

                 zij = zijl*th + zijr*(1.-th)             !- kpd - Jump in value
                 zid = zidl*th + zidr*(1.-th)             !- kpd - Jump in gradient. Equal to 0 here.

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1z=0 and rho2z=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Mixture density (Smeared)
                 rho1z(i,j,k+1) = rho1z(i,j,k+1)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE  
                 rho2z(i,j,k+1) = rho2z(i,j,k+1)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   - zij/aa/dz**2 - zid*th*(rho1/rho2)/aa/dz
                 w(i,j,k+1) = w(i,j,k+1) + zij/aa/dz**2 - zid*(1.-th)*(rho2/rho2)/aa/dz


                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigz(i,j,k+1) = - zij/aa/dz

                 icrv(i,j,k)   = 1
                 icrv(i,j,k+1) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell to front
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i,j,k+1).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k)*(1.-th) + crv(i,j,k+1)*th

                 zijl = xit*crv(i,j,k)
                 zijr = xit*crv(i,j,k+1)
                 zidl = 0.
                 zidr = 0.

                 zij = zijl*(1.-th) + zijr*th
                 zid = zidl*(1.-th) + zidr*th

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1z=0 and rho2z=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Density IS SMEARED HERE
                 rho1z(i,j,k+1) = rho1z(i,j,k+1)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2z(i,j,k+1) = rho2z(i,j,k+1)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + zij/aa/dz**2 + zid*(1.-th)*(rho2/rho2)/aa/dz
                 w(i,j,k+1) = w(i,j,k+1) - zij/aa/dz**2 + zid*th*(rho1/rho2)/aa/dz


                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigz(i,j,k+1) = zij/aa/dz

                 icrv(i,j,k) = 1
                 icrv(i,j,k+1) = 1

              end if

              !--------------------------------------------------------------
              !--------------------------------------------------------------

           end do
        end do
        end do

        crv = crv*icrv

      end subroutine mph_KPDcurvature3DC

!=========================================================================
!=========================================================================
!=========================================================================

