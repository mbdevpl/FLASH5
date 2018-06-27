subroutine sm_xdmfwrite(filebase,Nele,Nnodes,bodyType,Nodes, Elements)

#include "Flash.h"
#include "SolidMechanics.h"
  use Driver_interface, Only : Driver_abortFlash
  implicit none
  
  ! Argument list
  character(len=100),intent (IN) :: filebase
  integer, intent(IN) :: Nele, Nnodes, bodyType
  real, intent(IN) :: Nodes(Nnodes*NDIM);
  integer,intent(IN) :: Elements(NDIM,Nele)
  
  character(len=100) :: filename
  
end subroutine sm_xdmfwrite
