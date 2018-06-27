subroutine el10_lengths(XYZ,Qie,lengths)

  USE quadrule_C1_GL5
  USE sm_element_interface, only: el10_ShapeFunc

  implicit none    
#include "Flash.h"
#include "SolidMechanics.h"

  integer, parameter :: nen_e = NINE_NODES

  ! IO Variables
  real, dimension(nen_e,NDIM), intent(in)  :: XYZ, Qie
  real, dimension(2),          intent(out) :: lengths

  ! Internal variables
  integer, parameter :: nee = NDIM*nen_e
  real, external :: dnrm2, ddot
  integer :: i1, i2
  real    :: xi, eta
  real, dimension(nen_e):: NN, Ns, Nr
  real, dimension(NDIM) :: dr

  ! init the lengths
  lengths = 0.

  ! Compute length at (xi,eta) = (s,0)
  do i1 = 1,quad_nt
        xi = quad_xi(i1)
        eta = 0.
        
        call el10_ShapeFunc(NN,Nr,Ns,xi,eta)
        
        ! Compute the absolute position dr/dxi = Nr*(X + q)
        do i2 = 1,NDIM
           dr(i2) = ddot( nen_e, Nr, 1, XYZ(1:nen_e,i2)+Qie(1:nen_e,i2), 1)
        end do

        lengths(1) = lengths(1) + dnrm2(NDIM, dr, 1 )*quad_w(i1)

  end do

  ! Compute length at (xi,eta) = (0,s)
  do i1 = 1,quad_nt
        xi = 0.
        eta = quad_xi(i1)
        
        call el10_ShapeFunc(NN,Nr,Ns,xi,eta)
        
        ! Compute the absolute position dr/deta = Ns*(X + q)
        do i2 = 1,NDIM
           dr(i2) = ddot( nen_e, Ns, 1, XYZ(1:nen_e,i2)+Qie(1:nen_e,i2), 1)
        end do

        lengths(2) = lengths(2) + dnrm2(NDIM, dr, 1 )*quad_w(i1)

  end do
  

end subroutine el10_lengths
