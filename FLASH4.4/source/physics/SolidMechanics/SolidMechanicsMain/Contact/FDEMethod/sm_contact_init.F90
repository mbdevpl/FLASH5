  

subroutine sm_contact_init()
#include "Flash.h"

  use SolidMechanics_data, Only: sm_bodyInfo,  sm_NumBodies
  use sm_Contact_data, Only : sm_ContactInfo
  
  implicit none
  integer :: nnp,i


  if (sm_NumBodies>1) then
     allocate(sm_ContactInfo(sm_Numbodies))
     do i=1,sm_NumBodies
        nnp=sm_bodyInfo(i)%nnp
        allocate(sm_ContactInfo(i)%intersectedNodes(sm_NumBodies,nnp))
        allocate(sm_ContactInfo(i)%intersectForce(sm_NumBodies,NDIM))
        allocate(sm_ContactInfo(i)%intersect_count(sm_NumBodies))
        allocate(sm_ContactInfo(i)%intersect_vol(sm_NumBodies))
        sm_ContactInfo(i)%intersectedNodes(:,:) = 0
        sm_ContactInfo(i)%intersectForce(:,:)   = 0.0 
        sm_ContactInfo(i)%intersect_count(:) = 0
        sm_ContactInfo(i)%intersect_vol(:) = 0.0
     end do

  end if
  
  
end subroutine sm_contact_init

