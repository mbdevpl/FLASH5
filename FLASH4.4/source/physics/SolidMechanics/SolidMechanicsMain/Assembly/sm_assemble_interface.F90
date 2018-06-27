Module sm_assemble_interface

  implicit none

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

!
! Top Level:
!
  interface 
     subroutine sm_assemble_mass(ibd, iopt)
       implicit none
       integer, intent(IN) :: ibd
       integer, optional, intent(in) :: iopt
     end subroutine sm_assemble_mass
  end interface
  
  interface 
     subroutine sm_assemble_stiff(ibd, iopt)
       implicit none
       integer, intent(IN) :: ibd
       integer, optional, intent(in) :: iopt
     end subroutine sm_assemble_stiff
  end interface
  
  interface 
     subroutine sm_assemble_IntForce(ibd, iopt)
       implicit none
       integer, intent(IN) :: ibd
       integer, optional, intent(in) :: iopt
     end subroutine sm_assemble_IntForce
  end interface

  interface
     subroutine sm_assemble_ExtForce(ibd, time)
       implicit none
       integer, intent(in)    :: ibd! body number
       real, intent(in)       :: time
     end subroutine sm_assemble_ExtForce
  end interface

!
! 3D Flexible:
!
  interface 
     subroutine sm_assemble_mass_3DFlexible(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_mass_3DFlexible
  end interface
  
  interface 
     subroutine sm_assemble_stiff_3DFlexible(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_stiff_3DFlexible
  end interface
  
  interface 
     subroutine sm_assemble_IntForce_3DFlexible(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_IntForce_3DFlexible
  end interface

  interface
     subroutine sm_assemble_ExtForce_3DFlexible(ibd, time)
       implicit none
       integer, intent(in)    :: ibd! body number
       real, intent(in)       :: time
     end subroutine sm_assemble_ExtForce_3DFlexible
  end interface

!
! Rigid:
!
  interface 
     subroutine sm_assemble_mass_rigid(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_mass_rigid
  end interface
  
  interface 
     subroutine sm_assemble_stiff_rigid(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_stiff_rigid
  end interface
  
  interface 
     subroutine sm_assemble_IntForce_rigid(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_IntForce_rigid
  end interface

  interface
     subroutine sm_assemble_ExtForce_rigid(ibd, time)
       implicit none
       integer, intent(in)    :: ibd! body number
       real, intent(in)       :: time
     end subroutine sm_assemble_ExtForce_rigid
  end interface

!
! RBC
!
  interface 
     subroutine sm_assemble_mass_rbc(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_mass_rbc
  end interface
  
  interface 
     subroutine sm_assemble_stiff_rbc(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_stiff_rbc
  end interface
  
  interface 
     subroutine sm_assemble_IntForce_rbc(ibd, flag)
       implicit none
       integer, intent(IN) :: ibd, flag
     end subroutine sm_assemble_IntForce_rbc
  end interface

  interface
     subroutine sm_assemble_ExtForce_rbc(ibd, time)
       implicit none
       integer, intent(in)    :: ibd! body number
       real, intent(in)       :: time
     end subroutine sm_assemble_ExtForce_rbc
  end interface

  interface
     subroutine sm_assemble_COM(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_assemble_COM
  end interface

   interface
     subroutine sm_assemble_COM_3DFlexible(ibd, flag)
       implicit none
       integer, intent(in) :: ibd, flag
     end subroutine sm_assemble_COM_3DFlexible
  end interface

 interface
     subroutine sm_assemble_COM_rigid(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_assemble_COM_rigid
  end interface

 interface
     subroutine sm_assemble_COM_rbc(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_assemble_COM_rbc
  end interface
 
 interface
     subroutine sm_assemble_FluidVolInt_rigid()
       implicit none
     end subroutine sm_assemble_FluidVolInt_rigid
  end interface
 
  interface
     subroutine sm_assemble_FluidVolIntegrals()
       implicit none
     end subroutine sm_assemble_FluidVolIntegrals
  end interface



end Module sm_assemble_interface
