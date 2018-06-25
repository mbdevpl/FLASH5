
      subroutine mph_KPDlsRedistance(s,lsit,dx,dy,ix1,ix2,jy1,jy2)

        use IncompNS_data, ONLY : ins_cfl

        implicit none

#include "Flash.h"

        !- kpd - Imported variables

        integer, intent(in) :: lsit,ix1,ix2,jy1,jy2
        real, dimension(:,:,:), intent(inout):: s
        real, intent(in) :: dx,dy

        !- kpd - Local variables

        real :: so(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
                   sgn(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
                   eps,t,agf,dtL,err1,err2, &
                   ap,an,bp,bn,cp,cn,dp,dn, &
                   sxl,sxr,syl,syr,sm   

        integer :: i,j,k,n,incrm,m,itr

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
           s = 0.

           !- kpd - Loop through interior nodes
           do j = jy1,jy2
              do i = ix1,ix2

                 sxl = so(i-1,j,k)
                 sxr = so(i+1,j,k)
                 syl = so(i,j-1,k)
                 syr = so(i,j+1,k)

                 sm = so(i,j,k)

                 ap = max((sm-sxl),0.)/dx
                 an = min((sm-sxl),0.)/dx

                 bp = max((sxr-sm),0.)/dx
                 bn = min((sxr-sm),0.)/dx

                 cp = max((sm-syl),0.)/dy
                 cn = min((sm-syl),0.)/dy

                 dp = max((syr-sm),0.)/dy
                 dn = min((syr-sm),0.)/dy


                 if(so(i,j,k).gt.0.) then
                    agf = sqrt(max(ap**2,bn**2) + max(cp**2,dn**2))
                 elseif(so(i,j,k).lt.0.) then
                    agf = sqrt(max(an**2,bp**2) + max(cn**2,dp**2))
                 else
                    agf = 0.
                 end if

                 !solve Level Set distance equation -----------------------------

                 s(i,j,k) = so(i,j,k) + dtL*sgn(i,j,k)*(1. - agf)

                 err1 = max(err1,abs(1.-agf))

              end do
           end do

           err2 = sqrt(sum((s - so)**2)/dble(ix2-1)/dble(jy2-1))

!========================================================================
!========================================================================

           t = t + dtL

           if(itr.lt.lsit) then
              !write(6,'(i8,f10.5,2e20.10)') itr,t,err1
              print*,"LS Redistance Errors: ",itr,t,err1
              !m = m + incrm
           end if

        end do


      end subroutine mph_KPDlsRedistance 



