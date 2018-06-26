! Rigid assemble Center of Mass Function
! File:  
! Author: Marcos

#include "Flash.h"
#include "constants.h"

subroutine sm_assemble_COM_rigid(ibd)

  use SolidMechanics_data, ONLY: sm_BodyInfo

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd! body number
   
  ! Local variables:
  integer :: inod
  real :: xi,yi
#if NDIM == MDIM
  real :: zi
#endif

  ! Center of mass node:
  inod = sm_BodyInfo(ibd) % Borigin_node

  ! X - Y COM:
  xi= sm_BodyInfo(ibd)%x(inod) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(IAXIS,inod))
  yi= sm_BodyInfo(ibd)%y(inod) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(JAXIS,inod))

  sm_BodyInfo(ibd)%COM(IAXIS) = xi
  sm_BodyInfo(ibd)%COM(JAXIS) = yi  

  ! Z - COM:
#if NDIM == MDIM
  zi= sm_BodyInfo(ibd)%z(inod) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(KAXIS,inod))
  sm_BodyInfo(ibd)%COM(KAXIS) = zi
#endif

  return
    
end subroutine sm_assemble_COM_rigid
