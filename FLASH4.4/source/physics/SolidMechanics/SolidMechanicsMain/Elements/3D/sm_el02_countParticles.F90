!
! File:   sm_el10_countParticles.F90
! Author: tim
!

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_el02_countParticles(body, e, Dmin, nXi, nEta, ptelem, flag)
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
  integer, parameter :: nen_e = THREE_NODES
  real, dimension(nen_e,NDIM) :: XYZ, Qie, pt
  real :: Lseg, mXi1, mXi2
  real, external :: dnrm2

  ! get grid indo
  call get_Nodal_XYZ_ws(e, nen_e, XYZ, body)

  ! Get current Displacements
  select case( flag ) 
  case( SM_IOPT_QI )
     call get_Nodal_UVW_qi_ws(e, nen_e, Qie, body)
  case default
     ! This is Qn (not qi)
     call get_Nodal_UVW_qn_ws(e, nen_e, Qie, body)  
  end select 

  ! Position of triangles nodes:
  pt = XYZ + Qie

#ifdef TRIANG_MAPPING_A  /* Regular Paving */

  ! Compute nXi, X2-X1
  Lseg = dnrm2(NDIM, pt(2,1:NDIM) - pt(1,1:NDIM), 1)
  nXi = floor(Lseg/Dmin)
  !nXi = ceiling(Lseg/Dmin)

  ! Compute nEta, X3-X1
  Lseg = dnrm2(NDIM, pt(3,1:NDIM) - pt(1,1:NDIM), 1)
  nEta = floor(Lseg/Dmin)
  !nEta = ceiling(Lseg/Dmin)
  

  ! Get max for both
  nXi = max(nXi,nEta)
  if(nXi .le. 1 ) nXi = 1  ! The Infamous !!!
  nEta = nXi

  ! compute ptelem
  ptelem = (nXi**2 + nXi)/2

#else /* Double Loop Paving */

  ! Compute nXi, X2-X1
  Lseg = dnrm2(NDIM, pt(2,1:NDIM) - pt(1,1:NDIM), 1)
  nXi  = ceiling(2.*Lseg/(3.*Dmin)) ! Conservative estimate

  if(nXi .le. 1) nXi = 1
  nEta = nXi

  ! Compute ptelem  
  mXi1 = nXi*(nXi+1)/2
  mXi2 = (nXi-1)*nXi/2

  ptelem = mXi1 + mXi2

#endif


!  ! Compute nXi, X2-X1
!  pt = XYZ + Qie
!  Lseg = dnrm2(NDIM, pt(2,1:NDIM) - pt(1,1:NDIM), 1)
!  nXi = ceiling(Lseg/Dmin)
  
!  ! Compute nEta, X3-X1
!  Lseg = dnrm2(NDIM, pt(3,1:NDIM) - pt(1,1:NDIM), 1)
!  nEta = ceiling(Lseg/Dmin)

!  ! Get max for both
!  nXi = max(nXi,nEta)
!  if(nXi .eq. 1 ) nXi = 2  ! The Infamous !!!
!  nEta = nXi

!  ! compute ptelem
!  ptelem = (nXi**2 + nXi)/2 

  return

end subroutine sm_el02_countParticles
