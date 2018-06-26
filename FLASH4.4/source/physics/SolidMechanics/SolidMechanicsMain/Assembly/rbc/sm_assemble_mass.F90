!     
! File:   Assembly.F
! Author: tim
! Description: assembles mass matrix

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_assemble_mass(ibd)

  use SolidMechanics_data, only :  sm_nen, sm_BodyInfo, sm_structure
  use sm_element_interface, only: el05_mass, el12_mass
  use sm_Misc_interface, only: get_Nodal_XYZ

  use Driver_interface, only : Driver_abortFlash

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd ! body number
    
  ! Define internal variables
  type(sm_structure), pointer :: body
  integer :: nee, e, loop1, loop2, idx, nen_e
  real, allocatable, dimension(:,:) :: Me, XYZ
    
  body => sm_BodyInfo(ibd)
    
  ! Initialize M
  body % M   = 0.
  body % Mqv = 0.
    
  ! Loop over all the elements
  do e = 1,body%nel
        
     nen_e = sm_nen(body%eltype(e))
     nee = NDIM*nen_e
     allocate(Me(nee, nee))
        
     ! Get XYZ for the element
     allocate(XYZ(nen_e,3))
     call get_Nodal_XYZ(e, nen_e, XYZ, body)
        
     select case( body%eltype(e) )
     case( EIGHT_NODE_HEXAHEDRON )
        ! 8-Node hexahedron
        call el05_mass(Me, XYZ, body%MatDensity(e))
        
     case( TWSEVEN_NODE_HEXAHEDRON )
        ! 27-Node hexahedron
        call el12_mass(Me, XYZ, body%MatDensity(e) )
                
     case default
        call Driver_abortFlash("Error in sm_Assemble_stiff: unknown element type.")
     end select
    
     ! add into global mass matrix the same csc format as the stiffness matrix
     do loop1 = 1,nee
        
        ! check row
        if( body%LM(loop1,e) <= body%neq ) then
           do loop2 = 1,nee
              idx = body%LM_cs(loop1, loop2, e)
              
              ! Check Column:
              if( body%LM(loop2,e) <= body%neq ) then
                 body%M(idx) = body%M(idx) + Me(loop1,loop2)
              else
                 body%Mqv(idx) = body%Mqv(idx) + Me(loop1,loop2)
              endif
              
           enddo
        endif
        
     enddo
     
     deallocate(Me,XYZ)
     
  enddo

  return  
end subroutine sm_assemble_mass
