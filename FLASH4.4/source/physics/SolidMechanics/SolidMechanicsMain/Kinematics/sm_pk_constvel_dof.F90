!
!
!!!!
subroutine sm_pk_constvel_dof(time,maxrestparams,paramcoord,vc,vcd,vcdd)

#include "constants.h"
  implicit none
  integer, intent(in) :: maxrestparams
  real, intent(in)    :: time, paramcoord(maxrestparams)
  real, intent(out)   :: vc, vcd, vcdd
  
  real :: xo, vel

  ! Parameters: Are s.t. x(t) = xo + vel*t 
  xo = paramcoord(1)
  vel= paramcoord(2) 

  vc  = xo + vel*time 
  vcd =           vel 
  vcdd=           0.0 

  return

end subroutine sm_pk_constvel_dof
