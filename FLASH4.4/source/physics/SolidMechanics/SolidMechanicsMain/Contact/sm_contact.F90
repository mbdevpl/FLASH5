!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Contact/sm_contact
!!
!! NAME
!!  sm_Contact
!!
!! SYNOPSIS
!!
!! sm_Contact()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"

subroutine  sm_contact(ib,restart_local)

  use SolidMechanics_data, only : sm_meshMe

  implicit none

  ! Argument list
  integer, intent(IN) :: ib
  logical, intent(IN) :: restart_local
  
  if ((ib .eq. 1) .and. (sm_meshMe .eq. MASTER_PE)) write(*,*) 'THIS SIMULATION DOES NOT USE THE CONTACT UNIT..'

  return  
   
end subroutine sm_contact
 
