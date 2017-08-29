#include "Flash.h"

subroutine Grid_copyF4DataToMultiFabs(gds, phi, nodetype, reverse)

#ifdef FLASH_GRID_ANYAMREX
  use amrex_multifab_module
  implicit none
  type(amrex_multifab),OPTIONAL,intent(INOUT) :: phi(:)
#else
  type(*),OPTIONAL :: phi
#endif
  integer,intent(IN),OPTIONAL :: gds
  integer,intent(IN),OPTIONAL :: nodetype
  logical,intent(IN),OPTIONAL :: reverse

  ! Stub implementation, call is ignored!

end subroutine Grid_copyF4DataToMultiFabs
