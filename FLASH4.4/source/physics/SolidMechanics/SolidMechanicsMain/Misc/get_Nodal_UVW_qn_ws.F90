subroutine get_Nodal_UVW_qn_ws(e, nen_e, Qe, body)
  use SolidMechanics_data, only: sm_structure
  implicit none
#include "Flash.h"
#include "constants.h"    
  integer, intent(in) :: e, nen_e
  real, intent(out), dimension(nen_e,NDIM) :: Qe
  type(sm_structure) :: body
    
  integer :: a, idx1, idx2
    
  ! get qe = [ue, ve, we]
  do a = 1,nen_e
        
     idx1 = body%ws_IEN(a,e)
        
     idx2 = body%ID(IAXIS,idx1)
     Qe(a,IAXIS) = body%qn(idx2)
        
     idx2 = body%ID(JAXIS,idx1)
     Qe(a,JAXIS) = body%qn(idx2)

#if NDIM == MDIM        
     idx2 = body%ID(KAXIS,idx1)
     Qe(a,KAXIS) = body%qn(idx2)
#endif        
  end do
    
  return

end subroutine get_Nodal_UVW_qn_ws
