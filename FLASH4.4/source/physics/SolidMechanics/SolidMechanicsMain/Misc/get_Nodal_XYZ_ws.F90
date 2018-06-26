subroutine get_Nodal_XYZ_ws(e, nen_e, Xe, body)
  use SolidMechanics_data, only: sm_structure
  implicit none
#include "Flash.h"  
#include "constants.h"
  integer, intent(in) :: e, nen_e
  real, intent(out), dimension(nen_e,NDIM) :: Xe
  type(sm_structure), intent(in) :: body
   
  integer :: a, idx
    
  ! get Xe = [xe, ye, ze]
  do a = 1,nen_e
     idx = body%ws_IEN(a,e)  !on the wet surface
     Xe(a,IAXIS) = body%x(idx)
     Xe(a,JAXIS) = body%y(idx)
#if NDIM == MDIM
     Xe(a,KAXIS) = body%z(idx)
#endif
  enddo

  return
    
end subroutine get_Nodal_XYZ_ws
