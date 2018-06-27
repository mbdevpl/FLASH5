!
! File:   element_mass_matrices.F95
! Author: tim
!
! Created on February 2, 2012, 3:41 PM
!
#include "Flash.h"
#include "SolidMechanics.h"

subroutine el05_mass(me, XYZ, MatDensity)
    ! Compute element mass matrix for an 8 node hexahedron
    USE quadrule_C3_GL4
    USE sm_element_interface, only: el05_ShapeFunc, ShapeFunc_Expand

    implicit none

    integer, parameter :: nen_e = EIGHT_NODES
    integer, parameter :: nee = NDIM*nen_e
    real, dimension(nee,nee), intent(out) :: me
    real, dimension(nen_e,NDIM), intent(in) :: XYZ
    real, intent(in) :: MatDensity

    integer :: i1
    real  :: xi, eta, zeta, detJ, dV0
    real, dimension(nen_e)   :: NN, NNxi, NNeta, NNzeta
    real, dimension(NDIM,nee) :: N
    real, dimension(NDIM)     :: J1, J2, J3

    me = 0.

    ! Define Integration loop
    do i1 = 1,quad_nt
        xi = quad_xi(i1)
        eta = quad_eta(i1)
        zeta = quad_zeta(i1)

        ! Get Shape Functions, and local derivatives
        call el05_ShapeFunc(NN,NNxi,NNeta,NNzeta,xi,eta,zeta)

        ! Build [N]
        call ShapeFunc_Expand(NN, nen_e, N)

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

        ! Build me
        !me = me + rho(e)*(N'*N)*dV0
        call DGEMM('T','N', nee, nee, NDIM, MatDensity*dV0, N, & 
		            NDIM, N, NDIM, 1., me, nee)

    enddo

end
