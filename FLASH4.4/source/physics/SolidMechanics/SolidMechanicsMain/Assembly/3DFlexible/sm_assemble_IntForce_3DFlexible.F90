!     
! File:  
! Author: tim

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

subroutine sm_assemble_IntForce_3DFlexible(ibd, flag)

  use SolidMechanics_data, only: sm_nen,  sm_BodyInfo, sm_structure
  use sm_element_interface, only: el05_IntForce_Kirch, el05_IntForce_Biot, & 
                                  el12_IntForce_Kirch, el12_IntForce_Biot
  use sm_Misc_interface, only: get_Nodal_XYZ, get_Nodal_UVW_qn, get_Nodal_UVW_qi

  use Driver_interface, only : Driver_abortFlash

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd, flag ! body number
    
  ! Define internal variables
  type(sm_structure), pointer :: body
  integer :: nee, e, loop1, loop2, idx, nen_e
  real, allocatable, dimension(:,:) :: XYZe, Qie
  real, allocatable, dimension(:)   :: qes
    
  body => sm_BodyInfo(ibd)
  
  ! Initialize K and Q
  body % Qs  = 0.
    
  ! Loop over all the elements to calculate the area and volume terms
!$omp parallel do private(e, nee, nen_e, qes, XYZe, Qie, loop1, idx) shared(body)
  do e = 1, body%nel
     nen_e = sm_nen(body%eltype(e))  
     nee = NDIM*nen_e;
     allocate( qes(nee) )
     allocate(XYZe(nen_e,MDIM))
     allocate(Qie(nen_e,MDIM))
   
     ! get stuff:
     call get_Nodal_XYZ(e, nen_e, XYZe, body)

     select case( flag ) 
     case( SM_IOPT_QI )
        call get_Nodal_UVW_qi(e, nen_e, Qie, body)
     case default
        ! This is Qn (not qi)
        call get_Nodal_UVW_qn(e, nen_e, Qie, body)  
     end select
     
     select case( body%eltype(e) )

     case( EIGHT_NODE_HEXAHEDRON )
        ! 8-Node hexahedron
        select case( body%MatType(e) )
        case( MATERIAL_KIRCHHOFF )
           call el05_IntForce_Kirch( qes, XYZe, &
                         body%YoungsModulus(e), &
                         body%PoissonsRatio(e), Qie)
        case( MATERIAL_BIOT )
           call el05_IntForce_Biot( qes, XYZe, &
                        body%YoungsModulus(e), &
                        body%PoissonsRatio(e), Qie)
        end select
                
     case( TWSEVEN_NODE_HEXAHEDRON )
        ! 27-Node hexahedron
        select case( body%MatType(e) )
        case( MATERIAL_KIRCHHOFF )
           call el12_IntForce_Kirch( qes, XYZe, &
                body%YoungsModulus(e), &
                body%PoissonsRatio(e), Qie)
        case( MATERIAL_BIOT )
           call el12_IntForce_Biot( qes, XYZe, &
                body%YoungsModulus(e), &
                body%PoissonsRatio(e), Qie)
        end select
        
        
     case default
        call Driver_abortFlash("Error in sm_Assemble_IntForce_3DFlexible: unknown element type.")
     end select
     
     ! add into global internal force vector
     do loop1 = 1,nee
        idx = body%LM(loop1, e)
        ! check row:
        if( idx <= body%neq ) then
           body%Qs(idx) = body%Qs(idx) + qes(loop1)
        endif
     enddo        
     
     deallocate(qes, XYZe, Qie)      
  enddo
  
  return
    
end subroutine sm_assemble_IntForce_3DFlexible
