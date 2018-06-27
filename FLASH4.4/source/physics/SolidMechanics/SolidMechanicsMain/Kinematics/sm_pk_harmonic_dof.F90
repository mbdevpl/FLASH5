!
!
!!!!
subroutine sm_pk_harmonic_dof(time,maxrestparams,paramcoord,vc,vcd,vcdd)

#include "constants.h"
  implicit none
  integer, intent(in) :: maxrestparams
  real, intent(in)    :: time, paramcoord(maxrestparams)
  real, intent(out)   :: vc, vcd, vcdd
  
  real :: Ao,w,phase,fixed_coord


  ! Parameters: Are s.t. phi(t) = Ao*sin(w*t + phase) + fixed_coord
  Ao          =       paramcoord(1)
  w           = 2.*PI*paramcoord(2) ! Circular frequency
  phase       =       paramcoord(3)
  fixed_coord =       paramcoord(4)

  vc  =        Ao*sin(w*time + phase) + fixed_coord
  vcd =      w*Ao*cos(w*time + phase)
  vcdd=-(w**2)*Ao*sin(w*time + phase)

  return

end subroutine sm_pk_harmonic_dof
