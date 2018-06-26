!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/3D/sm_el10_ga_stiff
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!!  Linearize pressure:   p = p0 - h*dot( udd, nhat )
!!
!!  Qp = -\int N^T*p0*nhat*dS + \int N^T*h*dot(udd,nhat)*nhat*dS
!! 
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_el10_ga_stiff(body, e, particle, beta, dt, kf )
  use SolidMechanics_data,  only: sm_structure
  use sm_element_interface, only: el10_ShapeFunc, ShapeFunc_Expand
  use sm_Misc_interface,    only: sm_crossProd_mat
  use quadrule_C2_1
  implicit none

  ! constants
  integer, parameter :: nen_e = NINE_NODES
  integer, parameter :: nee   = NDIM*nen_e

  ! IO Variables
  type(sm_structure),           intent(in)  :: body
  integer,                      intent(in)  :: e  ! element number on the wet surface
  real, dimension(NPART_PROPS), intent(in)  :: particle
  real,                         intent(in)  :: beta, dt
  real, dimension(nee,nee),     intent(out) :: kf

  ! Internal Variables
  integer :: a, i,j,k, idx1, idx2, i1, i2, p, n_xi, n_eta, j1,j2, pL
  real, dimension(nen_e) :: NN, NNxi, NNeta
  real, dimension(NDIM, nee) :: N, Nxi, Neta, Dc
  real, dimension(nen_e, NDIM) :: pos, acc
  real, dimension(nee) :: qddn
  real :: xi, eta, r, s, xi1, xi2, eta1, eta2, dS, p0, nrmC, h
  real, dimension(NDIM) :: drdxi, drdeta, Nhat, c, udd
  real :: temp3x3(3,3), D2(3,nee), dndq(3,nee), temp1xnee(nee)
  real, external :: dnrm2, ddot

  !
  ! Get absoluate position and accl. of each node in the element
  !  p = a + nen_e*(i-1);
  do a = 1,nen_e
     idx1 = body%ws_IEN(a,e)
     pos(a,1) = body%x(idx1)
     pos(a,2) = body%y(idx1)
     pos(a,3) = body%z(idx1)
        
     do i = 1,3
        idx2     = body%ID(i,idx1)
        pos(a,i) = pos(a,i) + body%qn(idx2)
        acc(a,i) = body%qddn(idx2)
     end do
  enddo

  ! Get patch parameters
  n_xi  = body%ws_nXi(e)
  n_eta = body%ws_nEta(e)
  
  ! local patch number
  p = int( particle(PLOC_PART_PROP) ) 
  
  ! Map to local row/col info to define bounds
  j = ceiling( real(p)/real(n_xi) )
  k = p - (j-1)*n_xi
    
  ! compute bounds for integration
  xi1  = real(2*(k-1))/real(n_xi)  - 1.
  xi2  = real(2*(k-0))/real(n_xi)  - 1.
  eta1 = real(2*(j-1))/real(n_eta) - 1.
  eta2 = real(2*(j-0))/real(n_eta) - 1.

  ! get the pressure just outside the marker
  p0    = particle(PEXT_PART_PROP)

  ! get length h
  h = 1.2*particle(HL_PART_PROP)

  !
  ! Compute Integral
  !
  kf = 0.
  do i1 = 1,quad_nt
     r = quad_xi(i1)
     s = quad_eta(i1)
        
     ! Build (xi,eta)
     xi  = 0.5*(   (xi2 - xi1)*r +   (xi2 + xi1) )
     eta = 0.5*( (eta2 - eta1)*s + (eta2 + eta1) )
        
     call el10_ShapeFunc(NN,NNxi,NNeta,xi,eta)
     call ShapeFunc_Expand(NN,    nen_e, N)
     call ShapeFunc_Expand(NNxi,  nen_e, Nxi)
     call ShapeFunc_Expand(NNeta, nen_e, Neta)

     ! Compute surface Jacobian
     ! dx/dxi
     call DGEMV('T', nen_e, 3, 1., pos, nen_e, NNxi,   1, 0., drdxi, 1)
     ! dx/deta
     call DGEMV('T', nen_e, 3, 1., pos, nen_e, NNeta,  1, 0., drdeta, 1)
     c(1) = drdxi(2)*drdeta(3) - drdxi(3)*drdeta(2)
     c(2) = drdxi(3)*drdeta(1) - drdxi(1)*drdeta(3)
     c(3) = drdxi(1)*drdeta(2) - drdxi(2)*drdeta(1)
        
     ! compute dS
     nrmC = dnrm2(3,c,1)
     dS = 0.25*(xi2-xi1)*(eta2-eta1)*nrmC*quad_w(i1)

     ! Compute nhat
     nhat = c / nrmC

     ! compute dc/dq
     DC = matmul(sm_crossProd_mat(drdxi),Neta) - matmul(sm_crossProd_mat(drdeta), Nxi)

     ! compute the acceleration udd = [N]*{qdd}
     do i = 1,NDIM
        udd(i) = ddot(nen_e,NN,1,acc(:,i),1)
     end do

     ! compute $dnhat/dq = (1/|c|*I - (c.c)^-3/2*c*c^T)*A$
     temp3x3 = 0.
     do j1 = 1,3
        do j2 = 1,3
           temp3x3(j1,j2) = c(j1)*c(j2)
        end do
     end do
     temp3x3 = -1./nrmC**3 * temp3x3
     do j1 = 1,3
        temp3x3(j1,j1) = temp3x3(j1,j1) + 1./nrmC
     end do         
     dndq = matmul(temp3x3,DC)
     

     ! Compute $D2 = d/dq( \ddot{u}\cdot\hat{n} \vec{c} )$
     !         $   = c(n^T*N/(beta*dt^2) + udd^T*dn/dq)+ (udd.n)*DC $
     !
     ! compute n^T*N/(beta*dt^2)
     temp1xnee = 0.
     call daxpy(nen_e,nhat(1)/beta/dt**2, NN, 1, temp1xnee(1), 1)
     call daxpy(nen_e,nhat(2)/beta/dt**2, NN, 1, temp1xnee(nen_e+1),1)
     call daxpy(nen_e,nhat(3)/beta/dt**2, NN, 1, temp1xnee(2*nen_e+1),1)
     ! add it to udd*dndq
     temp1xnee = temp1xnee + matmul(udd,dndq)
     ! outer product with c
     do j1 = 1,3
        do j2 = 1,nee
           D2(j1,j2) = c(j1)*temp1xnee(j2)
        end do
     end do
     D2 = D2 + ddot(3,udd,1,nhat,1)*DC
     
     !Now compute the integral (this could be done faster since N is 2/3 zeros...maybe later.)
     kf = kf + matmul(transpose(N), -p0*DC + h*D2)*dS

  end do
     
end subroutine sm_el10_ga_stiff
