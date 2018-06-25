!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/sm_el10_mapParticles
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

subroutine sm_el10_mapParticles(body, e, ptelem,  &
                                xpos,ypos,zpos,   &
                                xvel,yvel,zvel,   &
                                xacc,yacc,zacc,   &
                                xnrm,ynrm,znrm,   &
                                areai, loc_num )
  use SolidMechanics_data, only: sm_structure
  use sm_element_interface, only: el10_ShapeFunc
  use quadrule_C2_1
  implicit none
  
  ! IO Variables
  type(sm_structure)   :: body     ! entire body structure
  integer, intent(in)  :: e        ! element number
  integer, intent(in)  :: ptelem
  real, dimension(ptelem) :: xpos, ypos, zpos, xvel, yvel, zvel, xacc, yacc, zacc, xnrm, ynrm, znrm, areai
  integer, dimension(ptelem) :: loc_num

  ! Internal variables
  integer, parameter :: nen_e = NINE_NODES
  integer :: a, i, j,k, idx1, idx2, i1, p, nxi, neta
  real, dimension(nen_e,NDIM) :: u, ud, udd
  real, dimension(nen_e) :: NN, NNxi, NNeta
  real :: xi, eta, r, s, xi1, xi2, eta1, eta2, area
  real, dimension(3) :: dxdxi, dxdeta, N
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
        ud(a,i)  = body%qdn(idx2)
        udd(a,i) = body%qddn(idx2)
     end do
  enddo

  nxi  = body%ws_nXi(e)
  neta = body%ws_nEta(e)

  !
  ! Loop over Particle Patches
  !
  do p = 1,ptelem

     ! set local patch number, cell centered on the natural square 
     ! 
     !  nEta   *    *       nXi*nEta
     !   *
     !   *                     *
     !   *
     !  nXi+1  *    *  
     !   1     2    3   *  *  nXi 
     !
     loc_num(p) = p

     ! Map to local row/col info to define bounds
     j = ceiling( real(p)/real(nxi) )
     k = p - (j-1)*nxi
    
     ! compute bounds for integration
     xi1  = real(2*(k-1))/real(nxi)  - 1.
     xi2  = real(2*(k-0))/real(nxi)  - 1.
     eta1 = real(2*(j-1))/real(neta) - 1.
     eta2 = real(2*(j-0))/real(neta) - 1.
     
     !
     ! Compute Area
     !
     area = 0
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
        
        ! compute area
        area = area + 0.25*(xi2-xi1)*(eta2-eta1)*dnrm2(3,N,1)*quad_w(i1)
     end do
     areai(p) = area
    
     ! Compute center point for Rpos, Rd, N calculations
     xi  = 0.5*( xi2 +  xi1)
     eta = 0.5*(eta2 + eta1) 
        
     ! Get shape functions
     call el10_ShapeFunc(NN,NNxi,NNeta,xi,eta)      
    
     ! Interp Position
     xpos(p) = ddot(nen_e, NN, 1, u(1:nen_e,1), 1)
     ypos(p) = ddot(nen_e, NN, 1, u(1:nen_e,2), 1)
     zpos(p) = ddot(nen_e, NN, 1, u(1:nen_e,3), 1)
     
     ! Interp Velocity
     xvel(p) = ddot(nen_e, NN, 1, ud(1:nen_e,1), 1)
     yvel(p) = ddot(nen_e, NN, 1, ud(1:nen_e,2), 1)
     zvel(p) = ddot(nen_e, NN, 1, ud(1:nen_e,3), 1)
    
     ! Interp Acceleration
     xacc(p) = ddot(nen_e, NN, 1, udd(1:nen_e,1), 1)
     yacc(p) = ddot(nen_e, NN, 1, udd(1:nen_e,2), 1)
     zacc(p) = ddot(nen_e, NN, 1, udd(1:nen_e,3), 1)
    
     !
     ! Compute normal vector
     !
     ! dx/dxi
     call DGEMV('T', nen_e, 3, 1., u, nen_e, NNxi,   1, 0., dxdxi, 1)
     ! dx/deta
     call DGEMV('T', nen_e, 3, 1., u, nen_e, NNeta,  1, 0., dxdeta, 1)
     N(1) = dxdxi(2)*dxdeta(3) - dxdxi(3)*dxdeta(2)
     N(2) = dxdxi(3)*dxdeta(1) - dxdxi(1)*dxdeta(3)
     N(3) = dxdxi(1)*dxdeta(2) - dxdxi(2)*dxdeta(1)
     N = N/dnrm2(3, N, 1) 
    
     xnrm(p) = N(1)
     ynrm(p) = N(2)
     znrm(p) = N(3)

  end do


  return

end subroutine sm_el10_mapParticles
