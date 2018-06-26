!
!
!
subroutine sm_el15_NBvars(dims,mdims,trmatrix,nTvar,Tpar,Tdpar,Tddpar,TNB,NwB_N,NaB_N) 

  use Driver_interface, only: Driver_abortFlash
  implicit none
#include "SolidMechanics.h"
  integer, intent(IN) :: dims,mdims,trmatrix,nTvar
  real, intent(IN) :: Tpar(nTvar),Tdpar(nTvar),Tddpar(nTvar)
  real, intent(OUT):: TNB(mdims,mdims),NwB_N(mdims,1),NaB_N(mdims,1)

  ! Local Variables

  if (dims .eq. mdims) then ! 3D Transformation

     select case (trmatrix)

     case(RB_IDENTITY)
        ! No rotation:
        TNB  =0.; TNB(1,1)=1.; TNB(2,2)=1.; TNB(3,3)=1.;
        NwB_N=0.;    NaB_N=0.;

     case(RB_EULER321)

        ! 3-2-1 Euler Angle Sequence



     case default
 
        call Driver_abortFlash("sm_el15_NBvars : Unknown Type of rigid body Transformation")

     end select

  else ! 2D one angle transformation


  endif

  return

end subroutine sm_el15_NBvars
