Module rbc_interface 

#include "constants.h"
#include "Flash.h"
  
  implicit none
  
  interface
     subroutine sm_rbc_Calc_parameters(ibd )
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_rbc_Calc_parameters
  end interface
  
end Module rbc_interface
