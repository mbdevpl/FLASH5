#include "ImBound.h"

subroutine ib_weightfunc(dsx,dsy,dsz,xp,yp,zp,x,y,z,wtype,wi,dwidx,dwidy,dwidz)

  implicit none
#include "Flash.h"
  integer, INTENT(IN) ::  wtype
  real, INTENT(IN) :: dsx,dsy,dsz,xp,yp,zp,x,y,z
  real, INTENT(OUT) :: wi,dwidx,dwidy,dwidz
      
  ! Local Variables
  real :: difx,dify,difz,adifx,adify,adifz,rix,riy,riz
  real :: drixdx,driydy,drizdz
  real :: wix,wiy,wiz,dwixdx,dwiydy,dwizdz

  real, parameter :: ep = 1.e-14

  difx = xp - x
  dify = yp - y
  adifx = abs(difx)
  adify = abs(dify)
  rix = adifx/dsx
  riy = adify/dsy

#ifdef IB_GET_DERIVS
  drixdx = 0.
  driydy = 0.

  if (adifx .gt. ep) then
     drixdx = sign(1.0,difx)/dsx
  endif

  if (adify .gt. ep) then
     driydy = sign(1.0,dify)/dsy
  endif
#endif

#if NDIM == 3
  difz  = zp - z
  adifz = abs(difz)
  riz   = adifz/dsz
#ifdef IB_GET_DERIVS
  drizdz = 0.
  if (adifz .gt. ep) then
     drizdz = sign(1.0,difz)/dsz
  endif
#endif
#else
  difz = 0.
  adifz= 0.
  riz  = 0. 
#endif

        
!!$  select case (wtype)
     
!!$  case (1) ! Cubic weight function distribution
        
     ! Wix and Wix,x:
     wix = 0.   
     dwixdx = 0.
     if (rix < 0.5) then
            
        wix = 2./3. - 4.*rix**2 + 4.*rix**3            
#ifdef IB_GET_DERIVS
        dwixdx = (-8.*rix + 12.*rix**2)*drixdx;
#endif            
     elseif (rix < 1.) then 
        
        wix = 4./3. - 4.*rix + 4.*rix**2 - 4./3.*rix**3
#ifdef IB_GET_DERIVS
        dwixdx = (-4. + 8.*rix - 4.*rix**2)*drixdx;
#endif            
     endif
            
     ! Wiy and Wiy,y:
     wiy = 0.
     dwiydy = 0.
     if (riy < 0.5) then
            
        wiy = 2./3. - 4.*riy**2 + 4.*riy**3;            
#ifdef IB_GET_DERIVS
        dwiydy = (-8.*riy + 12.*riy**2)*driydy;
#endif            
     elseif (riy < 1.) then
            
        wiy = 4./3. - 4.*riy + 4.*riy**2 - 4./3.*riy**3
#ifdef IB_GET_DERIVS
        dwiydy = (-4. + 8.*riy - 4.*riy**2)*driydy;
#endif
     endif
            
#if NDIM == 3

     ! Wiz and Wiz,z:
     wiz = 0.
     dwizdz = 0.
     if (riz < 0.5) then
        
        wiz = 2./3. - 4.*riz**2 + 4.*riz**3;
#ifdef IB_GET_DERIVS
        dwizdz = (-8.*riz + 12.*riz**2)*drizdz;
#endif
     elseif (riz < 1.) then

        wiz = 4./3. - 4.*riz + 4.*riz**2 - 4./3.*riz**3
#ifdef IB_GET_DERIVS
        dwizdz = (-4. + 8.*riz - 4.*riz**2)*drizdz;
#endif
     endif

#else
     wiz = 1.
#ifdef IB_GET_DERIVS
     dwizdz = 0.
#endif
#endif

     ! wi,dwidx,dwidy,dwidz:
     wi    = wix*wiy*wiz
#ifdef IB_GET_DERIVS
     dwidx = dwixdx*wiy*wiz
     dwidy = wix*dwiydy*wiz
     dwidz = wix*wiy*dwizdz
#else
     dwidx = 0.
     dwidy = 0.
     dwidz = 0.
#endif        

!!$  end select

  return

End Subroutine ib_weightfunc

