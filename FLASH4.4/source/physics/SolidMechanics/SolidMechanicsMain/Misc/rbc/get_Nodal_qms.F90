subroutine get_Nodal_qms(e, nen_e, Xe, body)
#include "Flash.h"
#include "SolidMechanics.h"
  use SolidMechanics_data, only: sm_structure
  implicit none
  
  integer, intent(in) :: e, nen_e
  real, intent(out), dimension(MAXNODERBC,NDIM) :: Xe
  type(sm_structure), intent(in) :: body
  
  integer :: a, idx1,idx2
  
 
  ! get Xe = [xe, ye, ze]
  do a = 1,nen_e
     idx1 = body%IEN(a,e)

     
     idx2 = body%ID(1,idx1)
     Xe(a,1) = body%qms(idx2,1)

     idx2 = body%ID(2,idx1)
     Xe(a,2) = body%qms(idx2,1)
    
#if NDIM==3
     idx2 = body%ID(3,idx1)
     Xe(a,3) = body%qms(idx2,1)

#endif
 
  enddo

end subroutine get_Nodal_qms
