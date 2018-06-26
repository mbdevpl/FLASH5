!
! File:   sm_el10_countParticles.F90
! Author: tim
!

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_el10_countParticles(body, e, Dmin, nXi, nEta, ptelem, flag)
  use SolidMechanics_data, only: sm_structure
  use sm_Misc_interface, only: get_Nodal_XYZ_ws, get_Nodal_UVW_qn_ws, get_Nodal_UVW_qi_ws
  use sm_element_interface, only: el10_lengths
  implicit none
  
  ! IO Variables
  type(sm_structure)   :: body     ! entire body structure
  integer, intent(in)  :: e        ! element number
  real, intent(in)     :: Dmin     ! Eulerian Distance
  integer, intent(out) :: nXi, nEta, ptelem ! number of points in Xi, Eta, and total
  integer, intent(in)  :: flag

  ! Internal variables
  integer, parameter :: nen_e = NINE_NODES
  real, dimension(nen_e,NDIM) :: XYZ, Qie
  real, dimension(2) :: lengths

  ! get grid info
  call get_Nodal_XYZ_ws(e, nen_e, XYZ, body)

  ! Get current Displacements
  select case( flag ) 
  case( SM_IOPT_QI )
     call get_Nodal_UVW_qi_ws(e, nen_e, Qie, body)  
  case default
     ! This is Qn (not qi)
     call get_Nodal_UVW_qn_ws(e, nen_e, Qie, body)  
  end select 

  ! compute lengths of the elements
  call el10_lengths(XYZ,Qie,lengths)

  ! Compute nXi
  nXi = ceiling(lengths(1)/Dmin)
  
  ! Compute nEta
  nEta = ceiling(lengths(2)/Dmin)

  ! compute ptelem
  ptelem = nXi*nEta

  return

end subroutine sm_el10_countParticles
