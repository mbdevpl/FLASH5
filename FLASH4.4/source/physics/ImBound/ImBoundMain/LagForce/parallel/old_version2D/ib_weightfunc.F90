
subroutine ib_weightfunc(dsx,dsy,xp,yp,x,y,wtype,wi,dwidx,dwidy)

  implicit none
  integer, INTENT(IN) ::  wtype
  real, INTENT(IN) :: dsx,dsy,xp,yp,x,y
  real, INTENT(OUT) :: wi,dwidx,dwidy
      
  ! Local Variables
  real :: difx,dify,adifx,adify,rix,riy
  real :: drixdx,driydy
  real :: wix,wiy,dwixdx,dwiydy
  real, parameter :: ep = 1.e-14

  difx = xp - x
  dify = yp - y

  adifx = abs(difx)
  adify = abs(dify)

  rix = adifx/dsx
  riy = adify/dsy

  drixdx = 0.
  driydy = 0.

  if (adifx > ep) then
     drixdx = sign(1.0,difx)/dsx
  endif

  if (adify > ep) then
     driydy = sign(1.0,dify)/dsy
  endif
        
  select case (wtype)
     
  case (1)
        
     ! Wix and Wix,x:
     if (rix < 0.5) then
            
        wix = 2./3. - 4.*rix**2 + 4.*rix**3            
        dwixdx = (-8.*rix + 12.*rix**2)*drixdx;
            
     elseif (rix < 1.) then 
        
        wix = 4./3. - 4.*rix + 4.*rix**2 - 4./3.*rix**3
        dwixdx = (-4. + 8.*rix - 4.*rix**2)*drixdx;
            
     else
        wix = 0.
        dwixdx = 0.
     endif
            
     ! Wiy and Wiy,y:
     if (riy < 0.5) then
            
        wiy = 2./3. - 4.*riy**2 + 4.*riy**3;            
        dwiydy = (-8.*riy + 12.*riy**2)*driydy;
            
     elseif (riy < 1.) then
            
        wiy = 4./3. - 4.*riy + 4.*riy**2 - 4./3.*riy**3
        dwiydy = (-4. + 8.*riy - 4.*riy**2)*driydy;
            
     else
        wiy = 0.
        dwiydy = 0.
     endif
            
     ! wi,dwidx,dwidy:
     wi = wix*wiy
     dwidx = dwixdx*wiy
     dwidy = wix*dwiydy
        
  end select

  return

End Subroutine ib_weightfunc

