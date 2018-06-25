1subroutine ib_extrapEpoints(xb,yb,sb,           &
                ielem,jelem,phile,zL,VarO,      &
                nx1,ny1,nz1,del,coord,bsize)


  use ImBound_data , only : ib_stencil

  implicit none
#include "Flash.h"
#include "constants.h"
  
!!! Argument list
  integer, INTENT(IN) :: nx1,ny1,nz1
  real, INTENT(IN) :: xb,yb,sb
  integer, INTENT(IN) :: ielem(ib_stencil),jelem(ib_stencil)
  real, INTENT(INOUT) :: phile(ib_stencil)
  real, INTENT(IN) :: zL
  real, INTENT(INOUT) :: VarO(nx1,ny1,nz1)
  real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)

  ! Local Variables:
  integer :: i
  integer :: ipos, jpos
  real :: dxloc, dyloc, hl, factor
  real :: dsp
  
  
  dsp = sb;
  ! Local cell size:
  dxloc = del(IAXIS)
  dyloc = del(JAXIS)

  ! Blocks associated Lagrangean points:
  
  !factor = 2.*dsp(lpindex(ii))/(dxloc+dyloc); ! Square grids 
  hl = 0.5*(dxloc+dyloc)                       ! Rectangular grids
  factor = dsp*hl/(dxloc*dyloc)

  phile(:) = factor*phile(:);
    
  do i = 1,ib_stencil
        
     ipos = ielem(i);
     jpos = jelem(i);
        
     VarO(ipos,jpos,1) = VarO(ipos,jpos,1) +  &
                         phile(i)*zL;
        

  enddo

  return

End Subroutine ib_extrapEpoints
