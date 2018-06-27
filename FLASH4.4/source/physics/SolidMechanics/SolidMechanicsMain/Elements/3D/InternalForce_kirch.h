!! IO: (qs,XYZ,YoungsModulus,PoissonRatio,Qie)

#if SM_ELEM_BUILDTYPE == EIGHT_NODE_HEXAHEDRON
    USE quadrule_C3_GL4
    USE sm_element_interface, only: el05_ShapeFunc, ShapeFunc_Expand
#elif SM_ELEM_BUILDTYPE == TWSEVEN_NODE_HEXAHEDRON
    USE quadrule_C3_GL5
    USE sm_element_interface, only: el12_ShapeFunc, ShapeFunc_Expand
#endif  

    implicit none    
#include "Flash.h"
#if SM_ELEM_BUILDTYPE == EIGHT_NODE_HEXAHEDRON
    integer, parameter :: nen_e = 8
#elif SM_ELEM_BUILDTYPE == TWSEVEN_NODE_HEXAHEDRON
    integer, parameter :: nen_e = 27
#endif  
	
    integer, parameter :: nee = NDIM*nen_e
    real, external :: DDOT
    
    real, dimension(nee), intent(out)     :: Qs
	
    real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
    real, intent(in) :: YoungsModulus,PoissonRatio
    
    integer :: i1
    real  :: xi, eta, zeta, detJ, dV0, mu, lambda
    real, dimension(3,3) :: Jinv
    real, dimension(nen_e)   :: NN, NNxi, NNeta, NNzeta, NNx, NNy, NNz
    real, dimension(nen_e,nen_e) :: dXX, dYY, dZZ, dYZ, dXZ, dXY, dYZ2, dXZ2, dXY2, tempD
    real, dimension(nen_e,3) :: DN
    real, dimension(nee,6)   :: Bt
    real, dimension(6,6)     :: C
    real, dimension(6)       :: S, Estrain
    real, dimension(nen_e,3)   :: XYZQe, XYZQe2, tempA
    
    ! Init qs
    qs = 0.
    
    ! get Xe and qe:
    !call get_Nodal_XYZ(e, nen_e, XYZ)
    !call get_Nodal_UVW_qi(e, nen_e, Qie)  
    XYZQe  = XYZ+Qie
    XYZQe2 = XYZ+0.5*Qie
        
    ! Get the material-constant matrix [C].
    mu = YoungsModulus/(2.*(1.+PoissonRatio))
    lambda = YoungsModulus*PoissonRatio/(1.+PoissonRatio)/(1.-2.*PoissonRatio)
    C = transpose(reshape( (/ &
             mu+mu+lambda,  lambda       ,  lambda       , 0., 0., 0., &
              lambda       , mu+mu+lambda,  lambda       , 0., 0., 0., &
              lambda       ,  lambda       , mu+mu+lambda, 0., 0., 0., &
                0.       ,    0.       ,    0.       ,   mu, 0., 0., &
                0.       ,    0.       ,    0.       , 0.,   mu, 0., &
                0.       ,    0.       ,    0.       , 0., 0.,   mu  /), (/6,6/)  ))    

    ! loop over the integral
    do i1 = 1,quad_nt
        xi = quad_xi(i1)
        eta = quad_eta(i1)
        zeta = quad_zeta(i1)
    
        ! Get Shape Functions, and local derivatives
#if SM_ELEM_BUILDTYPE == EIGHT_NODE_HEXAHEDRON
        call el05_ShapeFunc(NN,NNxi,NNeta,NNzeta,xi,eta,zeta)
#elif SM_ELEM_BUILDTYPE == TWSEVEN_NODE_HEXAHEDRON
        call el12_ShapeFunc(NN,NNxi,NNeta,NNzeta,xi,eta,zeta)
#endif      
        ! Compute Jacobian
        !J = [Nxi*Xe, Neta*Xe, Nzeta*Xe]
        !detJ = det(J)
        call compute_Jacobian3d(nen_e, XYZ, NNxi, NNeta, NNzeta, Jinv, detJ)
    
        ! Derivative matrix DN
        !DN = [NNxi', NNeta', NNzeta']/J
        tempA(1:nen_e,1) = NNxi
        tempA(1:nen_e,2) = NNeta
        tempA(1:nen_e,3) = NNzeta
        call DGEMM('N','N', nen_e, 3, 3, 1./detJ, tempA, nen_e, Jinv, 3, 0., DN, nen_e)                
    
        ! Build Nx, Ny, Nz
        NNx = DN(1:nen_e,1)
        NNy = DN(1:nen_e,2)
        NNz = DN(1:nen_e,3)
        
        ! Build dXX, dYY, dZZ, dYZ, dXZ, dXY
        !dXX = NNx'*NNx
        !dYY = NNy'*NNy
        !dZZ = NNz'*NNz
        !dYZ = NNy'*NNz
        !dXZ = NNx'*NNz
        !dXY = NNx'*NNy
        dXX = 0.
        call DGER( nen_e, nen_e, 1., NNx, 1, NNx, 1, dXX, nen_e)
        dYY = 0.
        call DGER( nen_e, nen_e, 1., NNy, 1, NNy, 1, dYY, nen_e)
        dZZ = 0.
        call DGER( nen_e, nen_e, 1., NNz, 1, NNz, 1, dZZ, nen_e)
        dYZ = 0.
        call DGER( nen_e, nen_e, 1., NNy, 1, NNz, 1, dYZ, nen_e)
        dXZ = 0.
        call DGER( nen_e, nen_e, 1., NNx, 1, NNz, 1, dXZ, nen_e)
        dXY = 0.
        call DGER( nen_e, nen_e, 1., NNx, 1, NNy, 1, dXY, nen_e)
        dYZ2 = Transpose(dYZ) + dYZ
        dXZ2 = Transpose(dXZ) + dXZ
        dXY2 = Transpose(dXY) + dXY
                       
        ! Build B
        !B = [ (X + qe)'*DXX
        !    (X + qe)'*DYY
        !    (X + qe)'*DZZ
        !    (X + qe)'*(DYZ' + DYZ)
        !    (X + qe)'*(DXZ' + DXZ)
        !    (X + qe)'*(DXY' + DXY)]
        tempA = 0.
        ! Part 1: (X + qe)'*DXX
        call DSYMM('L', 'L', nen_e, 3, 1., dXX, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,1) = reshape(tempA, (/ nee /) )       
        ! Part 2: (X + qe)'*DYY
        call DSYMM('L', 'L', nen_e, 3, 1., dYY, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,2) = reshape(tempA, (/ nee /) )
        ! Part 3: (X + qe)'*DZZ
        call DSYMM('L', 'L', nen_e, 3, 1., dZZ, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,3) = reshape(tempA, (/ nee /) )
        ! Part 4: (X + qe)'*(DYZ' + DYZ)
        call DSYMM('L', 'L', nen_e, 3, 1., dYZ2, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,4) = reshape(tempA, (/ nee /) )
        ! Part 5: (X + qe)'*(DXZ' + DXZ)
        call DSYMM('L', 'L', nen_e, 3, 1., dXZ2, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,5) = reshape(tempA, (/ nee /) )
        ! Part 6: (X + qe)'*(DXY' + DXY)
        call DSYMM('L', 'L', nen_e, 3, 1., dXY2, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,6) = reshape(tempA, (/ nee /) )
            
        !Estrain = 
        !  [ 1/2*qe'*DXX*qe + X'*DXX*qe  = (1/2*qe+X)'*DXX*qe
        !    1/2*qe'*DYY*qe + X'*DYY*qe
        !    1/2*qe'*DZZ*qe + X'*DZZ*qe
        !    qe'*DYZ*qe + X'*(DYZ'+DYZ)*qe
        !    qe'*DXZ*qe + X'*(DXZ'+DXZ)*qe
        !    qe'*DXY*qe + X'*(DXY'+DXY)*qe]
        ! Estrain(1) = (1/2*qe+X)'*DXX*qe:
        call DSYMM('L', 'L', nen_e, 3, 1., dXX, nen_e, Qie, nen_e, 0., tempA, nen_e)
        Estrain(1) = DDOT( nee, tempA, 1, XYZQe2, 1) 
        ! Estrain(2) = (1/2*qe+X)'*DYY*qe:
        call DSYMM('L', 'L', nen_e, 3, 1., dYY, nen_e, Qie, nen_e, 0., tempA, nen_e)
        Estrain(2) = DDOT( nee, tempA, 1, XYZQe2, 1) 
        ! Estrain(3) = (1/2*qe+X)'*DZZ*qe:
        call DSYMM('L', 'L', nen_e, 3, 1., dZZ, nen_e, Qie, nen_e, 0., tempA, nen_e)
        Estrain(3) = DDOT( nee, tempA, 1, XYZQe2, 1) 
        !
        ! Now the shears:
        !
        ! Estrain(4) = qe'*DYZ*qe + X'*(DYZ'+DYZ)*qe:
        call DSYMM('L', 'L', nen_e, 3, 1., dYZ, nen_e, Qie, nen_e, 0., tempA, nen_e)
        call DSYMM('L', 'L', nen_e, 3, 1., dYZ2, nen_e, XYZ, nen_e, 1., tempA, nen_e)
        Estrain(4) = DDOT( nee, tempA, 1, Qie, 1) 
        ! Estrain(5) = qe'*DXZ*qe + X'*(DXZ'+DXZ)*qe:
        call DSYMM('L', 'L', nen_e, 3, 1., dXZ, nen_e, Qie, nen_e, 0., tempA, nen_e)
        call DSYMM('L', 'L', nen_e, 3, 1., dXZ2, nen_e, XYZ, nen_e, 1., tempA, nen_e)
        Estrain(5) = DDOT( nee, tempA, 1, Qie, 1) 
        ! Estrain(6) = qe'*DXY*qe + X'*(DXY'+DXY)*qe:
        call DSYMM('L', 'L', nen_e, 3, 1., dXY, nen_e, Qie, nen_e, 0., tempA, nen_e)
        call DSYMM('L', 'L', nen_e, 3, 1., dXY2, nen_e, XYZ, nen_e, 1., tempA, nen_e)
        Estrain(6) = DDOT( nee, tempA, 1, Qie, 1) 
                               
        !S = C*Estrain
        call DGEMV('N', 6, 6, 1., C, 6, Estrain, 1, 0., S, 1)
    
        ! compute dV0
        dV0 = detJ*quad_w(i1)

        !
        ! Build Qs
        !
        ! Qs = Qs + B'*S*dV0
        call DGEMV( 'N', nee, 6, dV0, Bt, nee, S, 1, 1., qs, 1)
    
    enddo    
    
