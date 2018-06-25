!
!
!!!!
subroutine sm_pk_fixed_dof(time,maxrestparams,paramcoord,vc,vcd,vcdd)


  implicit none
  integer, intent(in) :: maxrestparams
  real, intent(in)    :: time, paramcoord(maxrestparams)
  real, intent(out)   :: vc, vcd, vcdd
  
  real :: fixed_coord

  fixed_coord = paramcoord(1)

  vc  = fixed_coord
  vcd = 0.
  vcdd= 0.

  return

end subroutine sm_pk_fixed_dof
