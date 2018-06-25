!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/sm_el10_FluidForce
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

subroutine sm_el10_FluidForce(body, e, particle, Hsp_pres, Hsp_visc )
  use SolidMechanics_data, only: sm_structure
  use quadrule_C2_1
  implicit none

  ! constants
  integer, parameter :: nen_e = NINE_NODES
  integer, parameter :: nee   = NDIM*nen_e

  ! IO Variables
  type(sm_structure),  intent(in) :: body
  integer, intent(in) :: e  ! element number on the wet surface
  real, dimension(NPART_PROPS), intent(in)  :: particle
  real, dimension(nee), intent(out) :: Hsp_pres, Hsp_visc

  ! Internal Variables
  integer :: a, i, j,k, idx1, idx2, i1, p, nxi, neta
  real, dimension(nen_e,NDIM) :: u
  real, dimension(nen_e) :: NN, NNxi, NNeta
  real :: xi, eta, r, s, xi1, xi2, eta1, eta2, dS
  real, dimension(NDIM) :: dxdxi, dxdeta, N, nhat, f_pres, f_visc
  real, external :: ddot, dnrm2

  !
  ! Get absoluate position of each node in the element
  !
  do a = 1,nen_e
     idx1 = body%ws_IEN(a,e)
     u(a,1) = body%x(idx1)
     u(a,2) = body%y(idx1)
     u(a,3) = body%z(idx1)
        
     do i = 1,3
        idx2     = body%ID(i,idx1)
        u(a,i)   = u(a,i) + body%qn(idx2)
     end do
  enddo

  ! Get patch parameters
  nxi  = body%ws_nXi(e)
  neta = body%ws_nEta(e)
  
  ! local patch number
  p = int( particle(PLOC_PART_PROP) ) 
  
  ! Map to local row/col info to define bounds
  j = ceiling( real(p)/real(nxi) )
  k = p - (j-1)*nxi
    
  ! compute bounds for integration
  xi1  = real(2*(k-1))/real(nxi)  - 1.
  xi2  = real(2*(k-0))/real(nxi)  - 1.
  eta1 = real(2*(j-1))/real(neta) - 1.
  eta2 = real(2*(j-0))/real(neta) - 1.

  ! get the unit norm of the patch
  nhat = (/ particle(NMLX_PART_PROP), particle(NMLY_PART_PROP), particle(NMLZ_PART_PROP) /)

  ! get the force on the patch from the particle
  ! $ f_pres = - pres*\hat{n} $
  f_pres = - particle(PRES_PART_PROP)*nhat
  ! $ f_visc = \vec{f}_\nu $
  f_visc = (/ particle(FXVI_PART_PROP), particle(FYVI_PART_PROP), particle(FZVI_PART_PROP) /)
     
  !
  ! Compute Integral
  !
  Hsp_pres(1:nee) = 0.
  Hsp_visc(1:nee) = 0.
  do i1 = 1,quad_nt
     r = quad_xi(i1)
     s = quad_eta(i1)
        
     ! Build (xi,eta)
     xi  = 0.5*(   (xi2 - xi1)*r +   (xi2 + xi1) )
     eta = 0.5*( (eta2 - eta1)*s + (eta2 + eta1) )
        
     call el10_ShapeFunc(NN,NNxi,NNeta,xi,eta)

     ! Compute surface Jacobian
     ! dx/dxi
     call DGEMV('T', nen_e, 3, 1., u, nen_e, NNxi,   1, 0., dxdxi, 1)
     ! dx/deta
     call DGEMV('T', nen_e, 3, 1., u, nen_e, NNeta,  1, 0., dxdeta, 1)
     N(1) = dxdxi(2)*dxdeta(3) - dxdxi(3)*dxdeta(2)
     N(2) = dxdxi(3)*dxdeta(1) - dxdxi(1)*dxdeta(3)
     N(3) = dxdxi(1)*dxdeta(2) - dxdxi(2)*dxdeta(1)
        
     ! compute dS
     dS = 0.25*(xi2-xi1)*(eta2-eta1)*dnrm2(3,N,1)*quad_w(i1)
     
     ! Compute $\int N^T* f_{pres} dS$
     ! Hsp = Hsp + N^T*f_pres*dS
     call daxpy(nen_e, f_pres(1)*dS, NN, 1, Hsp_pres(1), 1)
     call daxpy(nen_e, f_pres(2)*dS, NN, 1, Hsp_pres(nen_e+1), 1)
     call daxpy(nen_e, f_pres(3)*dS, NN, 1, Hsp_pres(2*nen_e+1), 1)

     ! Compute $\int N^T* f_{visc} dS$
     ! Hsp = Hsp + N^T*f_{visc}*dS
     call daxpy(nen_e, f_visc(1)*dS, NN, 1, Hsp_visc(1), 1)
     call daxpy(nen_e, f_visc(2)*dS, NN, 1, Hsp_visc(nen_e+1), 1)
     call daxpy(nen_e, f_visc(3)*dS, NN, 1, Hsp_visc(2*nen_e+1), 1)

  end do
     
end subroutine sm_el10_FluidForce
