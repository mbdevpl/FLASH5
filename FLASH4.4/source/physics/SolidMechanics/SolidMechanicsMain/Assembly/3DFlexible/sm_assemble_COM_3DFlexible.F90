!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Assembly/3DFlexible/sm_assemble_COM_3DFlexible
!!
!! NAME
!!  
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!!  
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

subroutine sm_assemble_COM_3DFlexible(ibd, flag)

  use SolidMechanics_data, only: sm_nen,  sm_BodyInfo, sm_structure
  use sm_element_interface, only: sm_el12_COM
  use sm_Misc_interface, only: get_Nodal_XYZ, get_Nodal_UVW_qn, get_Nodal_UVW_qi

  use Driver_interface, only : Driver_abortFlash

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd , flag! body number
    
  ! Define internal variables
  type(sm_structure), pointer :: body
  integer :: nee, e, nen_e
  real, allocatable, dimension(:,:) :: XYZe, Qie
  real, dimension(NDIM)   :: COM_e, moment_b
  real :: mass_e, mass_b
    
  body => sm_BodyInfo(ibd)
  
  ! Initalize the COM
  body % COM = 0.
  mass_b = 0.
  moment_b = 0.
    
  ! Loop over all the elements to calculate the area and volume terms
!$omp parallel do private(e, nen_e, nee, XYZe, Qie, mass_e, COM_e) shared(moment_b, mass_b, body, flag)
  do e = 1, body%nel
     nen_e = sm_nen(body%eltype(e))  
     nee = NDIM*nen_e
     allocate( XYZe(nen_e,MDIM) , Qie(nen_e,MDIM) )
     mass_e = 0.
     COM_e = 0.

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

        case( TWSEVEN_NODE_HEXAHEDRON )
           ! 27-Node hexahedron
           call sm_el12_COM(COM_e, mass_e, XYZe, Qie, body%MatDensity(e) )
          
        case default
           call Driver_abortFlash("Error in sm_Assemble_COM_3DFlexible: unknown element type.")
     end select
     
     ! sum up the mass and the moment
     moment_b = moment_b + mass_e*COM_e
     mass_b = mass_b + mass_e
     
     deallocate( XYZe, Qie )      
  enddo

  ! Compute the body centriod
  body%COM = moment_b/mass_b
  
  return
    
end subroutine sm_assemble_COM_3DFlexible
