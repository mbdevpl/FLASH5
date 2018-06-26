!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/sm_el12_COM
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
#include "SolidMechanics.h"

subroutine sm_el12_COM(COM, mass, XYZ, Qie,MatDensity)
! Compute element center of mass for a quadratic 27 node hexahedron

  USE quadrule_C3_GL5
  USE sm_element_interface, only: el12_ShapeFunc
  implicit none

  ! IO Variables
  integer, parameter :: nen_e = 27
  integer, parameter :: nee = NDIM*nen_e
  real, dimension(NDIM), intent(out) :: COM
  real                 , intent(out) :: mass
  real, dimension(nen_e,NDIM), intent(in) :: XYZ, Qie
  real, intent(in) :: MatDensity
	
  ! Internal Variables
  integer :: i1, i
  real    :: xi, eta, zeta, detJ, dV0
  real, dimension(nen_e)       :: NN, NNxi, NNeta, NNzeta
  real, dimension(nen_e, NDIM) :: XYZpQi
  real, dimension(NDIM)        :: J1, J2, J3, r, moment
  real, external :: ddot

  COM = 0.
  mass = 0.
  moment = 0.

  XYZpQi = XYZ + Qie

  ! Define Integration loop
  do i1 = 1,quad_nt
     xi = quad_xi(i1)
     eta = quad_eta(i1)
     zeta = quad_zeta(i1)

     ! Get Shape Functions, and local derivatives
     call el12_ShapeFunc(NN,NNxi,NNeta,NNzeta,xi,eta,zeta)
     
     ! Compute Jacobian
     !J = [Nxi*Xe, Neta*Xe, Nzeta*Xe]
     !detJ = det(J)
     call DGEMV('T', nen_e, 3, 1., XYZ, nen_e, NNxi,   1, 0., J1, 1)
     call DGEMV('T', nen_e, 3, 1., XYZ, nen_e, NNeta,  1, 0., J2, 1)
     call DGEMV('T', nen_e, 3, 1., XYZ, nen_e, NNzeta, 1, 0., J3, 1)
        
     detJ = +J1(1)*(J2(2)*J3(3) - J3(2)*J2(3)) &
            -J2(1)*(J1(2)*J3(3) - J3(2)*J1(3)) &
            +J3(1)*(J1(2)*J2(3) - J2(2)*J1(3))

     ! compute dV0
     dV0 = detJ*quad_w(i1)

     ! Build mass
     mass = mass + MatDensity*dV0

     ! Absolute position
     do i = 1,NDIM
        r(i) = ddot(nen_e, NN, 1, XYZpQi(1:nen_e, i), 1)
     end do

     ! Compute the moment $\int \vec{r} \rho_0 ~dV_0$
     moment = moment + r*MatDensity*dV0

  end do

  ! Build Center of Mass
  COM = moment/mass
    
end subroutine sm_el12_COM
