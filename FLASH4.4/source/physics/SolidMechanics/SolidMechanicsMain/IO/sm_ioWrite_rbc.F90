
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!          This subroutine write the RBC positions in ASCII format 
!          as Finite element data at each output step
!          written by Hussein Ezzeldin July 2010
!          future enhancments write the data in Binary format
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
subroutine sm_ioWrite_rbc(ibd,flag) !Positions,Elements,NNodes,NELE,NRBC)
  
#include "Flash.h"
#include "SolidMechanics.h"

  implicit none 
      
  ! Argument list
  integer, intent(IN) ::ibd
  integer, intent(IN), optional :: flag 
end subroutine sm_ioWrite_rbc
