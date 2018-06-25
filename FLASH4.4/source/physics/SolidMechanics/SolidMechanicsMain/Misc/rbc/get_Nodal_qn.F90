subroutine get_Nodal_qn(e, nen_e, Xe, body)
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
!!$     if (nen_e==4) then 
!!$        write(*,*) 'The bending point ',a,' is ', idx1,body%ID(1:3,idx1)
!!$        
!!$     end if
     idx2 = body%ID(1,idx1)
     Xe(a,1) = body%qn(idx2)
     
     idx2 = body%ID(2,idx1)
     Xe(a,2) = body%qn(idx2)
     if (NDIM==3) then
        idx2 = body%ID(3,idx1)
        Xe(a,3) = body%qn(idx2)
     end if
  enddo
  !if (nen_e==4) stop
end subroutine get_Nodal_qn
