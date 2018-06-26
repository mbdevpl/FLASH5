
      subroutine mph_KPDlsRedistance(s,u,v,dx,dy,ix1,ix2,jy1,jy2,soo,lsDT, blockID,minCellDiag)

        use IncompNS_data, ONLY : ins_cfl

        implicit none

#include "Flash.h"

        !- kpd - Imported variables
        integer, intent(in) :: ix1,ix2,jy1,jy2,blockID
        real, dimension(:,:,:), intent(inout):: s
        real, intent(in) :: dx,dy,lsDT, minCellDiag
        real, dimension(:,:,:), intent(in):: u,v,soo

        !- kpd - Local variables
        !real :: so(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
        !           sgn(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
        !           eps,t,agf,dtL,err1,err2, &
        !           ap,an,bp,bn,cp,cn,dp,dn, &
        !           sxl,sxr,syl,syr,sm,dt   
        real ::    so(NXB+2*NGUARD,NYB+2*NGUARD,1), &
                   sgn(NXB+2*NGUARD,NYB+2*NGUARD,1), &
                   eps,t,agf,dtL,err1,err2, &
                   ap,an,bp,bn,cp,cn,dp,dn, &
                   sxl,sxr,syl,syr,sm,dt,Dij   

        real :: phiXP,phiXM,phiYP,phiYM,xP,xM,yP,yM,eX

        integer :: i,j,k,n,m,itr

        !- kpd - For 2-D simulations
        k=1

        !t = 0.0
        eps = 1E-14
        eX = dx/1000.0
        itr = 0

        !- kpd - Level Set Redistance Iterations (pseudo-time)
        !do while(itr.lt.lsit)

        err1 = 0.

        itr = itr + 1

        !- kpd- Set the sign using the ORIGINAL distance function
        do j = jy1,jy2
           do i = ix1,ix2
              sgn(i,j,k) = soo(i,j,k)/abs(soo(i,j,k)+eps)
           end do
        end do

        !- kpd - 'so' and 's' are cell centered variables
        so  = s

        !- kpd - Loop through interior nodes
        do j = jy1,jy2
           do i = ix1,ix2

            s(i,j,k) = 0.
 
            !- kpd - Setup distance function stencil.'so' and 's' 
            !           are cell centered variables
            sm  = so(i,j,k)

            sxl = so(i-1,j,k)
            sxr = so(i+1,j,k)
            syl = so(i,j-1,k)
            syr = so(i,j+1,k)

            if ( (so(i,j,k)*so(i-1,j,k) .LT. 0) .OR. &
                 (so(i,j,k)*so(i+1,j,k) .LT. 0) .OR. &
                 (so(i,j,k)*so(i,j-1,k) .LT. 0) .OR. &
                 (so(i,j,k)*so(i,j+1,k) .LT. 0) ) then

               !=================================================================
               ! phi X+
               if ( (soo(i+1,j,k)*soo(i+2,j,k)).GT.0. .AND. (soo(i+1,j,k)*soo(i,j,k)).GT.0.) then
                  phiXP = soo(i,j,k) 
                  xP = dx/2.
               else
                  phiXP = soo(i+1,j,k)
                  xP = dx/2.
               end if               
               !=================================================================
               ! phiX-
               if ( (soo(i-1,j,k)*soo(i-2,j,k)).GT.0. .AND. (soo(i-1,j,k)*soo(i,j,k)).GT.0.) then
                  phiXP = soo(i,j,k) 
                  xM = dx/2.
               else
                  phiXP = soo(i-1,j,k)
                  xM = dx/2.
               end if               
               !=================================================================
               ! phi Y+
               if ( (soo(i,j+1,k)*soo(i,j+2,k)).GT.0. .AND. (soo(i,j+1,k)*soo(i,j,k)).GT.0.) then
                  phiYP = soo(i,j,k) 
                  yP = dy/2.
               else
                  phiYP = soo(i,j+1,k)
                  yP = dy/2.
               end if               
               !=================================================================
               ! phi Y-
               if ( (soo(i,j-1,k)*soo(i,j-2,k)).GT.0. .AND. (soo(i,j-1,k)*soo(i,j,k)).GT.0.) then
                  phiYP = soo(i,j,k) 
                  yM = dy/2.
               else
                  phiYP = soo(i,j-1,k)
                  yM = dy/2.
               end if               
               !=================================================================

               Dij = soo(i,j,k) / &
                     ( ( (phiXP-phiXM)/(MAX((xP-xM),(dx/1000.0))) )**2.0+ &
                       ( (phiYP-phiYM)/(MAX((yP-yM),(dx/1000.0))) )**2.0 )**0.5

               !------------------------------------------------------
               !- kpd - Solve Level Set Re-distance equation ---------
               !------------------------------------------------------
               s(i,j,k) = so(i,j,k) - lsDT/dx*(sgn(i,j,k)*(ABS(so(i,j,k)))-Dij)
               !------------------------------------------------------

            else

               !- KPD - First Order Upwind... +++++++++++++++++++++++
               ap = max((sm-sxl),0.)/dx
               an = min((sm-sxl),0.)/dx

               bp = max((sxr-sm),0.)/dx
               bn = min((sxr-sm),0.)/dx

               cp = max((sm-syl),0.)/dy
               cn = min((sm-syl),0.)/dy

               dp = max((syr-sm),0.)/dy
               dn = min((syr-sm),0.)/dy
               !+++++++++++++++++++++++++++++++++++++++++++++++++++++

               if(so(i,j,k).gt.0.) then
                  agf = ( sqrt(max(ap**2,bn**2) + max(cp**2,dn**2)) ) -1.0
               elseif(so(i,j,k).lt.0.) then
                  agf = ( sqrt(max(an**2,bp**2) + max(cn**2,dp**2)) ) -1.0
               else
                  agf = 0.
               end if
              
               !------------------------------------------------------
               !- kpd - Solve Level Set Re-distance equation ---------
               !------------------------------------------------------
               s(i,j,k) = so(i,j,k) - lsDT*(sgn(i,j,k)*(agf))
               !------------------------------------------------------

            end if  !if interface within one cell

            if (SIGN(1.0,s(i,j,k)) .NE. SIGN(1.0,soo(i,j,k))) then
               print*,"WARNING: LS Dist Function Changed Signs - ",i,j,k
            end if

           end do !do j = jy1,jy2
        end do    !do i = iy1,iy2

      end subroutine mph_KPDlsRedistance 

!========================================================================
!========================================================================

