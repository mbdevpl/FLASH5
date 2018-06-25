subroutine get_Nodal_UVW_qdm(e, nen_e, qde, body)
#include "Flash.h"
#include "SolidMechanics.h"
  
  use SolidMechanics_data, only: sm_structure
  implicit none
  
  integer, intent(in) :: e, nen_e
  real, intent(out), dimension(MAXNODERBC,NDIM) :: qde
  type(sm_structure) :: body
  
  integer :: a, idx1, idx2
   
  do a = 1,nen_e
     
     idx1 = body%IEN(a,e)
     
     idx2 = body%ID(1,idx1)
     qde(a,1) = body%qdms(idx2,1)
     
     idx2 = body%ID(2,idx1)
     qde(a,2) = body%qdms(idx2,1)
     if (NDIM==3) then
        idx2 = body%ID(3,idx1)
        qde(a,3) = body%qdms(idx2,1)
     end if
  enddo
  
end subroutine get_Nodal_UVW_qdm
