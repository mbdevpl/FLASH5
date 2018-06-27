!     
! File:   Assembly.F
! Author: tim
!
! Created on December 15, 2011, 12:10 PM
!

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_assemble_stiff(ibd)

  !USE struc_array, only: NDIM, nen, nel, eltype, MatType, max_eltype, K, Kqv, Qs, LM_cs, LM, neq
  
  use SolidMechanics_data, only: sm_nen,  sm_BodyInfo, sm_structure
  use sm_element_interface, only: el05_stiff_Kirch, el05_stiff_Biot, & 
                                  el12_stiff_Kirch, el12_stiff_Biot
  use sm_Misc_interface, only: get_Nodal_XYZ, get_Nodal_UVW_qi

  use Driver_interface, only : Driver_abortFlash


  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd ! body number
    
  ! Define internal variables
  type(sm_structure), pointer :: body
  integer :: nee, e, loop1, loop2, idx, nen_e
  real, allocatable, dimension(:,:) :: ke, XYZ, Qie
  real, allocatable, dimension(:)   :: qes
    
  body => sm_BodyInfo(ibd)
    
  ! Initialize K and Q
  body % K   = 0.
  body % Kqv = 0.
  body % Qs  = 0.
    
  ! Loop over all the elements
  do e = 1, body%nel
        
     nen_e = sm_nen(body%eltype(e))
     nee = NDIM*nen_e
     allocate( ke(nee, nee), qes(nee) )
        
     ! get stuff:
     call get_Nodal_XYZ(e, nen_e, XYZ, body)
     allocate(XYZ(nen_e,3), Qie(nen_e,3))
     call get_Nodal_UVW_qi(e, nen_e, Qie, body)  
    
                
     select case( body%eltype(e) )
     case( EIGHT_NODE_HEXAHEDRON )
        ! 8-Node hexahedron
        select case( body%MatType(e) )
        case( MATERIAL_KIRCHHOFF )
           call el05_stiff_Kirch( ke, qes, XYZ, &
                         body%YoungsModulus(e), &
                         body%PoissonsRatio(e), Qie)
        case( MATERIAL_BIOT )
           call el05_stiff_Biot( ke, qes, XYZ, &
                        body%YoungsModulus(e), &
                        body%PoissonsRatio(e), Qie)
        end select
                
     case( TWSEVEN_NODE_HEXAHEDRON )
        ! 27-Node hexahedron
        select case( body%MatType(e) )
        case( MATERIAL_KIRCHHOFF )
           call el12_stiff_Kirch( ke, qes, XYZ, &
                         body%YoungsModulus(e), &
                         body%PoissonsRatio(e), Qie)
        case( MATERIAL_BIOT )
           call el12_stiff_Biot( ke, qes, XYZ, &
                        body%YoungsModulus(e), &
                        body%PoissonsRatio(e), Qie)
        end select
                
     case default
        call Driver_abortFlash("Error in sm_Assemble_stiff: unknown element type.")
     end select
    
     ! add into global stiffness matrix, using csc format
     do loop1 = 1,nee
            
        idx = body%LM(loop1, e)
            
        ! check row:
        if( idx <= body%neq ) then
                
           body%Qs(idx) = body%Qs(idx) + qes(loop1)
           
           do loop2 = 1,nee
                    
              idx = body%LM_cs(loop1, loop2, e)
              ! check column:
              if( body%LM(loop2, e) <= body%neq ) then
                 body%K(idx) = body%K(idx) + ke(loop1,loop2)
              else
                 body%Kqv(idx) = body%Kqv(idx) + ke(loop1,loop2)
              endif
                    
           enddo
                
        endif
            
     enddo
        
     deallocate(ke, qes, XYZ, Qie)
        
  enddo

  return
    
end subroutine sm_assemble_stiff
