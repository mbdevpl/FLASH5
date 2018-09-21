!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_networkTable
!!
!! NAME
!!  
!!  bn_networkTable
!!
!!
!! SYNOPSIS
!! 
!!  call bn_networkTable()
!!
!!  
!! DESCRIPTION
!!
!!  routine networkTable generates the raw reaction rates for routine aprox13
!!  by using table interpolation instead of analytical expressions.
!!  a cubic polynomial is hardwired for speed
!!
!!***

subroutine bn_networkTable()

#include "Flash.h"

   
  use Burn_dataEOS, ONLY: btemp, bden
  use Burn_data
  use bn_dataAprox19

  implicit none

  !..uses tables instead of analytical expressions to evaluate the 
  !..raw reaction rates. a cubic polynomial is hardwired for speed.

  integer          i,j,imax,iat,mp,ifirst
  parameter        (mp = 4)
  real             tlo,thi,tstp,den_sav,btemp_sav,     &
       &                 x,x1,x2,x3,x4,a,b,c,d,e,f,g,h,p,q,  &
       &                 alfa,beta,gama,delt
  data             ifirst/0/


  !..make the table
  real :: entropy, dst, dsd
  if (ifirst .eq. 0) then
     ifirst = 1

     !..set the log temperature loop limits
     !..use 120 points per decade
     imax = 481
     if (imax .gt. nrattab) stop 'imax too small in bn_networkTable'
     tlo  = 6.0e0
     thi  = 10.0e0
     tstp = (thi - tlo)/float(imax-1)

     !..save the input
     btemp_sav = btemp
     den_sav   = bden

     !..form the table
     bden = 1.0e0
     do i=1,imax
        btemp = tlo + float(i-1)*tstp
        btemp = 10.0e0**(btemp)
        call bn_networkRates
        ttab(i) = btemp
        do j=1,nrat
           rattab(j,i) = ratraw(j)
        enddo
     enddo

     !..restore the input
     bden  = den_sav
     btemp = btemp_sav
  end if


  !..normal execution starts here
  !..set the density dependence vector
  dtab(irpp)    = bden
  dtab(irheng)  = bden
  dtab(irhegn)  = 1.0e0
  dtab(irhng)   = bden
  dtab(irdgn)   = 1.0e0
  dtab(irdpg)   = bden
  dtab(irhegp)  = 1.0e0
  dtab(ir33)    = bden
  dtab(ir34)    = bden
  dtab(ircpg)   = bden
  dtab(irnpg)   = bden
  dtab(ifa)     = 1.0e0
  dtab(ifg)     = 1.0e0
  dtab(iropg)   = bden
  dtab(irnag)   = bden
  dtab(ircag)   = bden
  dtab(iroga)   = 1.0e0
  dtab(ir3a)    = bden*bden
  dtab(irg3a)   = 1.0e0
  dtab(ir1212)  = bden
  dtab(ir1216)  = bden
  dtab(ir1616)  = bden
  dtab(iroag)   = bden
  dtab(irnega)  = 1.0e0
  dtab(irneag)  = bden
  dtab(irmgga)  = 1.0e0
  dtab(irmgag)  = bden
  dtab(irsiga)  = 1.0e0
  dtab(irmgap)  = bden
  dtab(iralpa)  = bden
  dtab(iralpg)  = bden
  dtab(irsigp)  = 1.0e0
  dtab(irsiag)  = bden
  dtab(irsga)   = 1.0e0
  dtab(irppa)   = bden
  dtab(irsiap)  = bden
  dtab(irppg)   = bden
  dtab(irsgp)   = 1.0e0
  dtab(irsag)   = bden
  dtab(irarga)  = 1.0e0
  dtab(irsap)   = bden
  dtab(irclpa)  = bden
  dtab(irclpg)  = bden
  dtab(irargp)  = 1.0e0
  dtab(irarag)  = bden
  dtab(ircaga)  = 1.0e0
  dtab(irarap)  = bden
  dtab(irkpa)   = bden
  dtab(irkpg)   = bden
  dtab(ircagp)  = 1.0e0
  dtab(ircaag)  = bden
  dtab(irtiga)  = 1.0e0
  dtab(ircaap)  = bden
  dtab(irscpa)  = bden
  dtab(irscpg)  = bden
  dtab(irtigp)  = 1.0e0
  dtab(irtiag)  = bden
  dtab(ircrga)  = 1.0e0
  dtab(irtiap)  = bden
  dtab(irvpa)   = bden
  dtab(irvpg)   = bden
  dtab(ircrgp)  = 1.0e0
  dtab(ircrag)  = bden
  dtab(irfega)  = 1.0e0
  dtab(ircrap)  = bden
  dtab(irmnpa)  = bden
  dtab(irmnpg)  = bden
  dtab(irfegp)  = 1.0e0
  dtab(irfeag)  = bden
  dtab(irniga)  = 1.0e0
  dtab(irfeap)  = bden
  dtab(ircopa)  = bden
  dtab(ircopg)  = bden
  dtab(irnigp)  = 1.0e0
  dtab(irfepg)  = bden
  dtab(ircogp)  = 1.0e0
  dtab(ir52ng)  = bden
  dtab(ir53gn)  = 1.0e0
  dtab(ir53ng)  = bden
  dtab(ir54gn)  = 1.0e0
  dtab(irpen)   = 1.0e0
  dtab(ispen)   = 1.0e0
  dtab(irnep)   = 1.0e0
  dtab(isnep)   = 1.0e0
  dtab(irn56ec) = 1.0e0
  dtab(isn56ec) = 1.0e0



  !..hash locate the temperature
  iat = int((log10(btemp) - tlo)/tstp) + 1
  iat = max(1,min(iat - mp/2 + 1,imax - mp + 1))

  !..setup the lagrange interpolation coefficients for a cubic
  x  = btemp
  x1 = ttab(iat)
  x2 = ttab(iat+1)
  x3 = ttab(iat+2)
  x4 = ttab(iat+3)
  a  = x - x1
  b  = x - x2
  c  = x - x3
  d  = x - x4
  e  = x1 - x2
  f  = x1 - x3
  g  = x1 - x4
  h  = x2 - x3
  p  = x2 - x4
  q  = x3 - x4
  alfa =  b*c*d/(e*f*g)
  beta = -a*c*d/(e*h*p)
  gama =  a*b*d/(f*h*q)
  delt = -a*b*c/(g*p*q)

  !..crank off the raw reaction rates
  do j=1,nrat
     ratraw(j) = (alfa*rattab(j,iat)    &
          &            + beta*rattab(j,iat+1)  &
          &            + gama*rattab(j,iat+2)  &
          &            + delt*rattab(j,iat+3)  &
          &              ) * dtab(j)
  enddo
  return
end subroutine bn_networkTable







