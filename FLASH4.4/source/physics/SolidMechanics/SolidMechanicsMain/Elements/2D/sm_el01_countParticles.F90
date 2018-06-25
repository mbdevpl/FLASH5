!
! File:   sm_el01_countParticles.F90
! 

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

subroutine sm_el01_countParticles(body, e, Dmin, nXi, nEta, ptelem, flag)
  use SolidMechanics_data, only: sm_structure
  use sm_Misc_interface, only: get_Nodal_XYZ_ws, get_Nodal_UVW_qn_ws, get_Nodal_UVW_qi_ws
  implicit none
  
  ! IO Variables
  type(sm_structure)   :: body     ! entire body structure
  integer, intent(in)  :: e        ! element number
  real, intent(in)     :: Dmin     ! Eulerian Distance
  integer, intent(out) :: nXi, nEta, ptelem ! number of points in Xi, Eta, and total
  integer, intent(in)  :: flag

  ! Internal variables
  integer, parameter :: nen_e = TWO_NODES
  real, dimension(nen_e,NDIM) :: XY, Qie, pt
  real :: Lseg, mXi1, mXi2
  real, external :: dnrm2

  ! get grid info
  call get_Nodal_XYZ_ws(e, nen_e, XY, body)

  ! Get current Displacements
  select case( flag ) 
  case( SM_IOPT_QI )
     call get_Nodal_UVW_qi_ws(e, nen_e, Qie, body)
  case default
     ! This is Qn (not qi)
     call get_Nodal_UVW_qn_ws(e, nen_e, Qie, body)  
  end select 

  ! Position of triangles nodes:
  pt = XY + Qie

  ! Compute nXi, X2-X1
  Lseg = dnrm2(NDIM, pt(2,1:NDIM) - pt(1,1:NDIM), 1)
  nXi = ceiling(Lseg/Dmin)

  ! nEta is not used on lines:
  nEta = CONSTANT_ONE

  ! compute ptelem
  ptelem = nXi

  return

end subroutine sm_el01_countParticles
