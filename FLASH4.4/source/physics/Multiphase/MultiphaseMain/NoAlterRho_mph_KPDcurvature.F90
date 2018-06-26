!=========================================================================
!=========================================================================
!=========================================================================


        subroutine mph_KPDcurvature2DAB(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,visc,vis1,vis2)

   
        implicit none

#include "Flash.h"

        !*************************************************************************

        !---------------------------------
        !- kpd - Data from routine call...
        !---------------------------------
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx, dy, rho1, rho2, xit, vis1, vis2 
        real, intent(out) :: crmx, crmn
        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,sigx,sigy,visc

        !--------------------------
        !- kpd - Local variables...
        !--------------------------
        integer :: i,j,k
        real :: rPhiXN,rPhiXE,rPhiXS,rPhiXW, &
                rPhiYN,rPhiYE,rPhiYS,rPhiYW, &
                rMagN,rMagE,rMagS,rMagW

        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr
        real, parameter :: eps = 1E-13

        !*************************************************************************

        !- kpd - For 2D runs 
        k=1

        !- kpd - Compute the curvature ---------------

        crv = 0.
        do j = jy1,jy2
           do i = ix1,ix2    
              !----------------------------------------------------
              !- kpd - 2 phi gradients per face method
              !----------------------------------------------------
              !        X - Location
              rPhiXE = 1./dx*(s(i+1,j,k)-s(i,j,k)  )
              rPhiXW = 1./dx*(s(i,j,k)  -s(i-1,j,k))
              rPhiXN = 1./4./dx * ( (s(i+1,j+1,k) - s(i-1,j+1,k)) &
                                  + (s(i+1,j,k)   - s(i-1,j,k)  ) )
              rPhiXS = 1./4./dx * ( (s(i+1,j,k)   - s(i-1,j,k)  ) &
                                  + (s(i+1,j-1,k) - s(i-1,j-1,k)) )
              !        Y - Location
              rPhiYN = 1./dy*(s(i,j+1,k)-s(i,j,k)  )
              rPhiYS = 1./dy*(s(i,j,k)  -s(i,j-1,k))
              rPhiYE = 1./4./dy * ( (s(i+1,j+1,k) - s(i+1,j-1,k)) &
                                  + (s(i,j+1,k)   - s(i,j-1,k)  ) )
              rPhiYW = 1./4./dy * ( (s(i,j+1,k)   - s(i,j-1,k)  ) &
                                  + (s(i-1,j+1,k) - s(i-1,j-1,k)) )
              !----------------------------------------------------

              !- kpd - Compute the magnitude of the gradient at each face
              rMagE = sqrt( rPhiXE**2. + rPhiYE**2. ) + eps
              rMagW = sqrt( rPhiXW**2. + rPhiYW**2. ) + eps
              rMagN = sqrt( rPhiXN**2. + rPhiYN**2. ) + eps
              rMagS = sqrt( rPhiXS**2. + rPhiYS**2. ) + eps


              crv(i,j,k) = 1./dx * (rPhiXE/rMagE - rPhiXW/rMagW) &
                         + 1./dy * (rPhiYN/rMagN - rPhiYS/rMagS)
              !----------------------------------------------------

!if (i.eq.4 .and. j.eq.4)then
!!print*,rPhiYN,rMagN,rPhiYS,rMagS
!!print*,s(i,j+1,k),s(i,j,k),s(i,j-1,k)
!print*,"Curv Val",i,j,k,crv(i,j,k)
!end if

           end do
        end do

        !*************************************************************************

        !----------------------------------
        !- kpd - Compute the phase function
        !----------------------------------
        k=1
        do j = jy1-1,jy2+1
           do i = ix1-1,ix2+1
              pf(i,j,k) = 0.

              if(s(i,j,k).ge.0.) then
                 pf(i,j,k) = 1.                       !- kpd - Set phase function on each side of interface
                 visc(i,j,k) = vis1/vis2               !- kpd - Set viscosity on each side of interface
                 !print*,"vis1",i,j,visc(i,j,k)
              else
                 visc(i,j,k) = vis2/vis2
                 !print*,"vis2",i,j,visc(i,j,k)
              end if
              !print*,"VISC",i,j,k,visc(i,j,k)
              !print*,"PFUN",i,j,pf(i,j,k),visc(i,j,k)
           end do
        end do

        !--------------------------------------------------------------
        !- kpd - These are FACE VALUED inverse densities for each phase
        !--------------------------------------------------------------

        !- kpd - density on x-face
        rho1x = 0.
        rho2x = 0.
        !- kpd - Loop through boundary and interior cell faces
        do j = jy1-1,jy2+1
           do i = ix1-1,ix2+1
              a1 = (pf(i-1,j,k) + pf(i,j,k)) / 2.                       
              a2 = pf(i-1,j,k)  /abs(pf(i-1,j,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)
              !rho1x(i,j,k) = a1*a2/rho1
              !rho2x(i,j,k) = (1. - a1*a2)/rho2
              rho1x(i,j,k) = a1*a2/(rho1/rho2)
              rho2x(i,j,k) = (1. - a1*a2)/(rho2/rho2)
           end do
        end do

        !- kpd - density on y-face
        rho1y = 0.
        rho2y = 0.
        !- kpd - Loop through boundary and interior cell faces
        do i = ix1-1,ix2+1
           do j = jy1-1,jy2+1
              a1 = (pf(i,j-1,k) + pf(i,j,k)) / 2.           
              a2 = pf(i,j-1,k)  /abs(pf(i,j-1,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)
              !rho1y(i,j,k) = a1*a2/rho1
              !rho2y(i,j,k) = (1. - a1*a2)/rho2
              rho1y(i,j,k) = a1*a2/(rho1/rho2)
              rho2y(i,j,k) = (1. - a1*a2)/(rho2/rho2)
           end do
        end do



      end subroutine mph_KPDcurvature2DAB

!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================

!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!- kpd - There is a break in the ...curvature routine to 
!        go back and do guard cell filling for curvature.
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------

!=========================================================================
!=========================================================================
!=========================================================================
!- kpd - From here on mph_curvature2DC:
!=========================================================================
!=========================================================================
!=========================================================================

        subroutine mph_KPDcurvature2DC(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,refLev)   

   
        implicit none

#include "Flash.h"

        !implicit real*8(a-h,o-z)
        !common/param/re,xit,rho1,rho2,g,sp,ubc,uout,nint

        integer, intent(in) :: ix1,ix2,jy1,jy2, refLev
        real, intent(in) :: dx, dy, rho1, rho2, xit 
        real, intent(out) :: crmx, crmn

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,sigx,sigy
        integer :: icrv(NXB+2*NGUARD,NYB+2*NGUARD,1)

        !- kpd - 
        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr
        integer :: i,j,k
        real, parameter :: eps = 1E-13

        integer :: iSmear

!--------------------------------------------
!----------------jump conditions ------------
!--------------------------------------------
!xij: jump in value
!xid: jump in gradient 
!l,r, values at pts left and right of the interface 
!crv: curvature
!cri: curvature at interface 
!left interface between i and i+1
!  phase 2 at i, phase 1 at i+1
!  theta = x_i+1 - x_I
!right interface between i and i+1
!  phase 1 at i, phase 2 at i+1
!  theta = x_I - x_i  
!--------------------------------------------
!--------------------------------------------

!print*,refLev

        iSmear  = 1

        crmx = -1E10
        crmn = 1E10
        sigx = 0.
        sigy = 0.
        w = 0.
        icrv = 0

        !- kpd - Need to loop through one guard cell on each side to set jumps 
        !           when they cross block boundaries
                    !do j = jy1,jy2
                    !   do i = ix1,ix2
        k=1
        do j = jy1-1,jy2
           do i = ix1-1,ix2
        !do j = jy1-1,jy2+1
        !   do i = ix1-1,ix2+1


              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i+1,j,k).eq.1.) then

                 !th = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i,j,k))) 

!                 if (refLev .eq. 1) then
!                    th = 0.5
!                 else
!                    th = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i,j,k)))
!                 end if
                 

                 cri = crv(i+1,j,k)*(1.-th) + crv(i,j,k)*th

                 xijl = xit*crv(i,j,k)                 !- kpd - sigma*K. Used for jump in pressure
                 xijr = xit*crv(i+1,j,k)               !- kpd - sigma*K. Used for jump in pressure
                 xidl = 0.                             !- kpd - Used for jump in gradient
                 xidr = 0.                             !- kpd - Used for jump in gradient

                 xij = xijl*th + xijr*(1.-th)          !- kpd - Jump in value
                 xid = xidl*th + xidr*(1.-th)          !- kpd - Jump in gradient. Equal to 0 here.
                 
                 if (iSmear  .eq. 1) then
!- kpd - SMEARED OUT JUMPS...
                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)           !- kpd - Mixture density (not inverse)
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa 

!                 if (refLev .eq. 1) then
!                    rho1x(i+1,j,k) = 0.0
!                    rho2x(i+1,j,k) = 0.5*(rho1+rho2)
!                 end if

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   - xij/aa/dx**2 - xid*th*(rho1/rho2)/aa/dx  
                 w(i+1,j,k) = w(i+1,j,k) + xij/aa/dx**2 - xid*(1.-th)*(rho2/rho2)/aa/dx
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = - xij/aa/dx           !- kpd - sigma*K/rho/dx 
         
                 else
!- kpd - NO SMEARING OUT JUMPS (SHARP)
                 w(i,j,k)   = w(i,j,k)   - (((rho1x(i,j,k)  +rho2x(i,j,k) )) *rho2)*xit*crv(i,j,k)  /1.0/dx**2   
                 w(i+1,j,k) = w(i+1,j,k) + (((rho1x(i+1,j,k)+rho2x(i+1,j,k)))*rho2)*xit*crv(i+1,j,k)/1.0/dx**2 
                 sigx(i+1,j,k) = - (((rho1x(i+1,j,k)+rho2x(i+1,j,k)))*rho2)*xit*crv(i+1,j,k)/1.0/dx
                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i+1,j,k).eq.0.) then

!                 if (refLev .eq. 1) then
!                    th = 0.5
!                 else
!                    th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
!                 end if

                 cri = crv(i,j,k)*(1.-th) + crv(i+1,j,k)*th

                 xijl = xit*crv(i,j,k)
                 xijr = xit*crv(i+1,j,k)
                 xidl = 0. 
                 xidr = 0.

                 xij = xijl*(1.-th) + xijr*th
                 xid = xidl*(1.-th) + xidr*th

                 if (iSmear  .eq. 1) then
!- kpd - SMEARING OUT JUMPS 
                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa   

!                 if (refLev .eq. 1) then
!                    rho1x(i+1,j,k) = 0.0
!                    rho2x(i+1,j,k) = 0.5*(rho1+rho2)
!                 end if

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + xij/aa/dx**2 + xid*(1.-th)*(rho2/rho2)/aa/dx
                 w(i+1,j,k) = w(i+1,j,k) - xij/aa/dx**2 + xid*th*(rho1/rho2)/aa/dx
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = xij/aa/dx

                 else
!- kpd - NO SMEARING OUT JUMPS (SHARP)
                 w(i,j,k)   = w(i,j,k)   + (((rho1x(i,j,k)  +rho2x(i,j,k)  ))*rho2)*xit*crv(i,j,k)  /1.0/dx**2   
                 w(i+1,j,k) = w(i+1,j,k) - (((rho1x(i+1,j,k)+rho2x(i+1,j,k)))*rho2)*xit*crv(i+1,j,k)/1.0/dx**2 
                 sigx(i+1,j,k) = (((rho1x(i+1,j,k)+rho2x(i+1,j,k)))*rho2)*xit*crv(i+1,j,k)/1.0/dx
                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i,j+1,k).eq.1.) then

!                 if (refLev .eq. 1) then
!                    th = 0.5
!                 else
!                    th = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j,k)))
!                 end if

                 cri = crv(i,j+1,k)*(1.-th) + crv(i,j,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0. 
                 yidr = 0.

                 yij = yijl*th + yijr*(1.-th)
                 yid = yidl*th + yidr*(1.-th)

                 if (iSmear  .eq. 1) then
!- kpd - SMEARING OUT JUMPS 
                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa 

!                 if (refLev .eq. 1) then
!                    rho1y(i+1,j,k) = 0.0
!                    rho2y(i+1,j,k) = 0.5*(rho1+rho2)
!                 end if

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k) - yij/aa/dy**2 - yid*th*(rho1/rho2)/aa/dy
                 w(i,j+1,k) = w(i,j+1,k)   + yij/aa/dy**2 - yid*(1.-th)*(rho2/rho2)/aa/dy 
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = - yij/aa/dy

                 else
!- kpd - NO SMEARING OUT JUMPS (SHARP)
                 w(i,j,k)   = w(i,j,k)   - (((rho1y(i,j,k)  +rho2y(i,j,k)  )) *rho2)*xit*crv(i,j,k)   /1.0/dy**2
                 w(i,j+1,k) = w(i,j+1,k) + (((rho1y(i,j+1,k)+rho2y(i,j+1,k))) *rho2)*xit*crv(i,j+1,k) /1.0/dy**2
                 sigy(i,j+1,k) = - (((rho1y(i,j+1,k)+rho2y(i,j+1,k)))*rho2) *xit*crv(i,j+1,k)/1.0/dy
                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i,j+1,k).eq.0.) then

!                 if (refLev .eq. 1) then
!                    th = 0.5
!                 else
!                    th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))) 
!                 end if
 
                 cri = crv(i,j,k)*(1.-th) + crv(i,j+1,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0. 
                 yidr = 0.

                 yij = yijl*(1.-th) + yijr*th
                 yid = yidl*(1.-th) + yidr*th

                 if (iSmear  .eq. 1) then
!- kpd - SMEARING OUT JUMPS 
                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa  

!                 if (refLev .eq. 1) then
!                    rho1y(i+1,j,k) = 0.0
!                    rho2y(i+1,j,k) = 0.5*(rho1+rho2)
!                 end if

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + yij/aa/dy**2 + yid*(1.-th)*(rho2/rho2)/aa/dy
                 w(i,j+1,k) = w(i,j+1,k) - yij/aa/dy**2 + yid*th*(rho1/rho2)/aa/dy  
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = yij/aa/dy

!!print*,"TH",i,j,th,crv(i,j,k),crv(i,j+1,k),w(i,j,k),w(i,j+1,k),sigy(i,j+1,k)

                 else
!- kpd - NO SMEARING OUT JUMPS (SHARP)
                 w(i,j,k)   = w(i,j,k)   + (((rho1y(i,j,k)  +rho2y(i,j,k)  )) *rho2)*xit*crv(i,j,k)   /1.0/dy**2
                 w(i,j+1,k) = w(i,j+1,k) - (((rho1y(i,j+1,k)+rho2y(i,j+1,k))) *rho2)*xit*crv(i,j+1,k) /1.0/dy**2
                 sigy(i,j+1,k) = (((rho1y(i,j+1,k)+rho2y(i,j+1,k)))*rho2) *xit*crv(i,j+1,k)/1.0/dy
                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !--------------------------------------------------------------

              !print*,"MPHsigp",i,j,w(i,j,k),pf(i,j,k),pf(i+1,j,k)
              !print*,"MPH",i,j,rho1x(i,j,k)+rho2x(i,j,k),rho1y(i,j,k)+rho2y(i,j,k)

              !print*,"CurvDens",i,j,rho1x(i,j,k)+rho2x(i,j,k),rho1y(i,j,k)+rho2y(i,j,k)

           end do
        end do

        !- kpd - This is done for post-processing to visualize where jumps were applied...
        crv = crv*icrv

      end subroutine mph_KPDcurvature2DC


