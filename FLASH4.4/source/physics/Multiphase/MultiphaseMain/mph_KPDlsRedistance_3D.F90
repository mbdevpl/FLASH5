
      subroutine mph_KPDlsRedistance_3D(s,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,soo,lsDT, minCellDiag)

        use IncompNS_data, ONLY : ins_cfl

        implicit none

#include "Flash.h"

        !- kpd - Imported variables

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, dimension(:,:,:), intent(inout):: s
        real, intent(in) :: dx,dy,dz,lsDT, minCellDiag
        real, dimension(:,:,:), intent(in):: u,v,w,soo

        !- kpd - Local variables

        !real :: so(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
        !           sgn(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
        real :: so(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
                   sgn(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
                   eps,t,agf,dtL,err1,err2, &
                   ap,an,bp,bn,cp,cn,dp,dn, &
                   sxl,sxr,syl,syr,sm, &
                   ep,en,fp,fn,szl,szr,Dij   

        integer :: i,j,k,n,incrm,m

        eps = 1E-14
        incrm = 1
        m = 1
        err1 = 0.

        !- kpd- Set the sign using the ORIGINAL distance function
        do k = kz1,kz2
           do j = jy1,jy2
              do i = ix1,ix2
                 sgn(i,j,k) = soo(i,j,k)/abs(soo(i,j,k)+eps)
              end do
           end do
        end do

        !- kpd - 'so' and 's' are cell centered variables
        so = s

        !- kpd - Loop through interior nodes
        !        'so' and 's' are cell centered variables
        do k = kz1,kz2
           do j = jy1,jy2
              do i = ix1,ix2

                 s(i,j,k) = 0.
 
                 !- kpd - If cell is near interface...
                 if ( (so(i,j,k)*so(i-1,j,k) .LT. 0) .OR. &
                      (so(i,j,k)*so(i+1,j,k) .LT. 0) .OR. &
                      (so(i,j,k)*so(i,j-1,k) .LT. 0) .OR. &
                      (so(i,j,k)*so(i,j+1,k) .LT. 0) .OR. &
                      (so(i,j,k)*so(i,j,k-1) .LT. 0) .OR. &
                      (so(i,j,k)*so(i,j,k+1) .LT. 0) ) then


                    !------------------------------------------------------
                    !- kpd - Solve Level Set Re-distance equation ---------
                    !------------------------------------------------------
                    s(i,j,k) = soo(i,j,k) 
                    !------------------------------------------------------


                 !- kpd - If cell is NOT near interface...
                 else

                   !- kpd - Setup distance function stencil.
                   sxl = so(i-1,j,k)
                   sxr = so(i+1,j,k)
                   syl = so(i,j-1,k)
                   syr = so(i,j+1,k)
                   szl = so(i,j,k-1)
                   szr = so(i,j,k+1)

                   sm = so(i,j,k)

                   ap = max((sm-sxl),0.)/dx
                   an = min((sm-sxl),0.)/dx

                   bp = max((sxr-sm),0.)/dx
                   bn = min((sxr-sm),0.)/dx

                   cp = max((sm-syl),0.)/dy
                   cn = min((sm-syl),0.)/dy

                   dp = max((syr-sm),0.)/dy
                   dn = min((syr-sm),0.)/dy

                   ep = max((sm-szl),0.)/dz
                   en = min((sm-szl),0.)/dz

                   fp = max((szr-sm),0.)/dz
                   fn = min((szr-sm),0.)/dz

                   !---------------------------------------------
                   !- kpd - Compute the magnitude of the gradient
                   !---------------------------------------------
                   if(soo(i,j,k).gt.0.) then
                      agf = sqrt( max(ap**2,bn**2) &
                                + max(cp**2,dn**2) &
                                + max(ep**2,fn**2) ) - 1.0
                   elseif(soo(i,j,k).lt.0.) then
                      agf = sqrt( max(an**2,bp**2) &
                                + max(cn**2,dp**2) &
                                + max(en**2,fp**2) ) - 1.0
                   else
                      agf = 0.
                   end if

                   !------------------------------------------------------
                   !- kpd - Solve Level Set Re-distance equation ---------
                   !------------------------------------------------------
                   s(i,j,k) = so(i,j,k) - lsDT*(sgn(i,j,k)*(agf))
                   !------------------------------------------------------

                 end if

                 if (SIGN(1.0,s(i,j,k)) .NE. SIGN(1.0,soo(i,j,k))) then
                    print*,"WARNING: LS Dist Function Changed Signs - ",i,j,k
                 end if

                 end do    !do i = ix1,ix2
              end do       !do j = jx1,jx2
           end do          !do k = kx1,kx2

           err2 = sqrt(sum((s - so)**2)/dble(ix2-1)/dble(jy2-1)/dble(kz2-1))

      end subroutine mph_KPDlsRedistance_3D 



