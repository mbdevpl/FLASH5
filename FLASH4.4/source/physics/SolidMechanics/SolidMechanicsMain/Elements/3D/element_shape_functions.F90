!     
! File:   element_shape_functions.F90
! Author: tim
!
! Shape function for finite element meshes generated using GMSH node ordering
! 
#include "SolidMechanics.h"

subroutine el01_ShapeFunc(NN,Nr,r)
! Compute shape functions for 2 node line
    implicit none
    real, intent(in)  :: r
    integer, parameter  :: n = TWO_NODES
    real, intent(out) :: NN(n), Nr(n)
    
    NN = (/ 0.5*(1.-r), 0.5*(1.+r)  /)
    Nr = (/ -0.5, 0.5 /)
        
  end subroutine el01_ShapeFunc

subroutine el02_ShapeFunc(NN,Nr,Ns,r,s)
! Compute shape functions for 3 node triangle
    implicit none
    real, intent(in)  :: r,s
    integer, parameter  :: n = THREE_NODES
    real :: u
    real, intent(out) :: NN(n), Nr(n), Ns(n)
    
    u = 1.d0 - r - s
    
    NN = (/ u, r, s /)
    Nr = (/ -1.d0, 1.d0, 0.d0 /)
    Ns = (/ -1.d0, 0.d0, 1.d0 /)
    
  end subroutine el02_ShapeFunc
    
subroutine el05_ShapeFunc(NN,Nr,Ns,Nt,r,s,t)
! Compute shape functions for a linear 8 node hexahedron
    implicit none
    real, intent(in)  :: r,s,t
    integer, parameter  :: n = EIGHT_NODES
    real, intent(out) :: NN(n), Nr(n), Ns(n), Nt(n)
    
    NN = (/-((-1. + r)*(-1. + s)*(-1. + t))/8.,      &
        ((1. + r)*(-1. + s)*(-1. + t))/8.,           &
        -((1. + r)*(1. + s)*(-1. + t))/8.,           &
        ((-1. + r)*(1. + s)*(-1. + t))/8.,           &
        ((-1. + r)*(-1. + s)*(1. + t))/8.,           &
        -((1. + r)*(-1. + s)*(1. + t))/8.,           &
        ((1. + r)*(1. + s)*(1. + t))/8.,             &
        -((-1. + r)*(1. + s)*(1. + t))/8. /)
    Nr = (/-((-1. + s)*(-1. + t))/8.,                &
        ((-1. + s)*(-1. + t))/8.,                    &
        -((1. + s)*(-1. + t))/8.,                    &
        ((1. + s)*(-1. + t))/8.,                     &
        ((-1. + s)*(1. + t))/8.,                     &
        -((-1. + s)*(1. + t))/8.,                    &
        ((1. + s)*(1. + t))/8.,                      &
        -((1. + s)*(1. + t))/8. /)
    Ns = (/-((-1. + r)*(-1. + t))/8.,                &
        ((1. + r)*(-1. + t))/8.,                     &
        -((1. + r)*(-1. + t))/8.,                    &
        ((-1. + r)*(-1. + t))/8.,                    &
        ((-1. + r)*(1. + t))/8.,                     &
        -((1. + r)*(1. + t))/8.,                     &
        ((1. + r)*(1. + t))/8.,                      &
        -((-1. + r)*(1. + t))/8. /)
    Nt = (/ -((-1. + r)*(-1. + s))/8.,               &
        ((1. + r)*(-1. + s))/8.,                     &
        -((1. + r)*(1. + s))/8.,                     &
        ((-1. + r)*(1. + s))/8.,                     &
        ((-1. + r)*(-1. + s))/8.,                    &
        -((1. + r)*(-1. + s))/8.,                    &
        ((1. + r)*(1. + s))/8.,                      &
        -((-1. + r)*(1. + s))/8. /)
    
end

subroutine el09_ShapeFunc(NN,Nr,Ns,r,s)
! Compute shape functions for 6 node triangle
    implicit none
    real, intent(in)  :: r,s
    integer, parameter  :: n = SIX_NODES
    real :: u
    real, intent(out) :: NN(n), Nr(n), Ns(n)
    
    u = 1. - r - s
    
    NN = (/ u*(-1. + 2.*u),      &
        r*(-1. + 2.*r),          &
        s*(-1. + 2.*s),          &
        4.*r*u, 4.*r*s, 4.*s*u/)
    Ns = (/ 1. - 4.*u, -1. + 4.*r, 0.d0 , -4.*r + 4.*u, 4.*s, -4.*s /)
    Nr = (/ 1. - 4.*u, 0.d0 , -1. + 4.*s, -4.*r, 4.*r, -4.*s + 4.*u /)
    
  end subroutine el09_ShapeFunc

subroutine el10_ShapeFunc(NN,Nr,Ns,r,s)
! Compute shape functions for 9 node quad

    implicit none
    real, intent(in)  :: r,s
    integer,parameter   :: n = NINE_NODES
    real, intent(out) :: NN(n), Nr(n), Ns(n)
    
    NN = (/ ((-1.0 + r)*r*(-1.0 + s)*s)/4.,       &
        (r*(1.0 + r)*(-1.0 + s)*s)/4.,            &
        (r*(1.0 + r)*s*(1.0 + s))/4.,             &
        ((-1.0 + r)*r*s*(1.0 + s))/4.,            &
        -((-1.0 + r)*(1.0 + r)*(-1.0 + s)*s)/2.,  &
        -(r*(1.0 + r)*(-1.0 + s)*(1.0 + s))/2.,   &
        -((-1.0 + r)*(1.0 + r)*s*(1.0 + s))/2.,   &
        -((-1.0 + r)*r*(-1.0 + s)*(1.0 + s))/2.,  &
        (-1.0 + r)*(1.0 + r)*(-1.0 + s)*(1.0 + s) /)
    Nr = (/ ((-1.0 + 2.0*r)*(-1.0 + s)*s)/4.,     &
        ((1.0 + 2.0*r)*(-1.0 + s)*s)/4.,          &
        ((1.0 + 2.0*r)*s*(1.0 + s))/4.,           &
        ((-1.0 + 2.0*r)*s*(1.0 + s))/4.,          &
        -(r*(-1.0 + s)*s),                        &
        -((1.0 + 2.0*r)*(-1.0 + s**2))/2.,        &
        -(r*s*(1.0 + s)),                         &
        -((-1.0 + 2.0*r)*(-1.0 + s**2))/2.,       &
        2.0*r*(-1.0 + s**2) /)
    Ns = (/ ((-1.0 + r)*r*(-1.0 + 2.0*s))/4.,     &
        (r*(1.0 + r)*(-1.0 + 2.0*s))/4.,          &
        (r*(1.0 + r)*(1.0 + 2.0*s))/4.,           &
        ((-1.0 + r)*r*(1.0 + 2.0*s))/4.,          &
        -((-1.0 + r**2)*(-1.0 + 2.0*s))/2.,       &
        -(r*(1.0 + r)*s),                         &
        -((-1.0 + r**2)*(1.0 + 2.0*s))/2.,        &
        -((-1.0 + r)*r*s),                        &
        2.0*(-1.0 + r**2)*s /)
        
  end subroutine el10_ShapeFunc

subroutine el11_ShapeFunc(NN,Nr,Ns,Nt,r,s,t)
! Compute shape functions for a 10-node second order tetrahedron 
!(4 nodes associated with the vertices and 6 with the edges)

    implicit none 
    real, intent(in)  :: r,s,t
    integer,parameter   :: n = TEN_NODES
    real              :: u
    real, intent(out) :: NN(n), Nr(n), Ns(n), Nt(n)
    
    u = 1. - r - s - t
    
    NN = (/ u*(-1. + 2.*u), r*(-1. + 2.*r), s*(-1. + 2.*s), t*(-1. + 2.*t), &
        4.*r*u, 4.*r*s, 4.*s*u, 4.*t*u, 4.*s*t, 4.*r*t /)
    Nr = (/ 1. - 4.*u, -1. + 4.*r, 0d0, 0d0, -4.*r + 4.*u, 4.*s, -4.*s, -4.*t, 0d0, 4.*t /)
    Ns = (/ 1. - 4.*u, 0d0, -1. + 4.*s, 0d0, -4.*r, 4.*r, -4.*s + 4.*u, -4.*t, 4.*t, 0d0 /)
    Nt = (/ 1. - 4.*u, 0d0, 0d0, -1. + 4.*t, -4.*r, 0d0, -4.*s, 4.*u - 4.*t, 4.*s, 4.*r /)
    
  end subroutine el11_ShapeFunc
    
    

subroutine el12_ShapeFunc(NN,Nr,Ns,Nt,r,s,t)
! Compute shape functions for a quadratic 27 node hexahedron

    implicit none 
    real, intent(in)  :: r,s,t
    integer,parameter   :: n = TWENTYSEVEN_NODES
    real, intent(out) :: NN(n), Nr(n), Ns(n), Nt(n)
        
    NN = (/ ((-1.0 + r)*r*(-1.0 + s)*s*(-1.0 + t)*t)/8.,            (r*(1.0 + r)*(-1.0 + s)*s*(-1.0 + t)*t)/8.,                  &
        (r*(1.0 + r)*s*(1.0 + s)*(-1.0 + t)*t)/8.,                  ((-1.0 + r)*r*s*(1.0 + s)*(-1.0 + t)*t)/8.,                  &
        ((-1.0 + r)*r*(-1.0 + s)*s*t*(1.0 + t))/8.,                 (r*(1.0 + r)*(-1.0 + s)*s*t*(1.0 + t))/8.,                   &
        (r*(1.0 + r)*s*(1.0 + s)*t*(1.0 + t))/8.,                   ((-1.0 + r)*r*s*(1.0 + s)*t*(1.0 + t))/8.,                   &
        -((-1.0 + r)*(1.0 + r)*(-1.0 + s)*s*(-1.0 + t)*t)/4.,       -((-1.0 + r)*r*(-1.0 + s)*(1.0 + s)*(-1.0 + t)*t)/4.,        &
        -((-1.0 + r)*r*(-1.0 + s)*s*(-1.0 + t)*(1.0 + t))/4.,        -(r*(1.0 + r)*(-1.0 + s)*(1.0 + s)*(-1.0 + t)*t)/4.,        &
        -(r*(1.0 + r)*(-1.0 + s)*s*(-1.0 + t)*(1.0 + t))/4.,         -((-1.0 + r)*(1.0 + r)*s*(1.0 + s)*(-1.0 + t)*t)/4.,        &
        -(r*(1.0 + r)*s*(1.0 + s)*(-1.0 + t)*(1.0 + t))/4.,          -((-1.0 + r)*r*s*(1.0 + s)*(-1.0 + t)*(1.0 + t))/4.,        &
        -((-1.0 + r)*(1.0 + r)*(-1.0 + s)*s*t*(1.0 + t))/4.,         -((-1.0 + r)*r*(-1.0 + s)*(1.0 + s)*t*(1.0 + t))/4.,        &
        -(r*(1.0 + r)*(-1.0 + s)*(1.0 + s)*t*(1.0 + t))/4.,          -((-1.0 + r)*(1.0 + r)*s*(1.0 + s)*t*(1.0 + t))/4.,         &
        ((-1.0 + r)*(1.0 + r)*(-1.0 + s)*(1.0 + s)*(-1.0 + t)*t)/2., ((-1.0 + r)*(1.0 + r)*(-1.0 + s)*s*(-1.0 + t)*(1.0 + t))/2.,&
        ((-1.0 + r)*r*(-1.0 + s)*(1.0 + s)*(-1.0 + t)*(1.0 + t))/2., (r*(1.0 + r)*(-1.0 + s)*(1.0 + s)*(-1.0 + t)*(1.0 + t))/2., &
        ((-1.0 + r)*(1.0 + r)*s*(1.0 + s)*(-1.0 + t)*(1.0 + t))/2.,  ((-1.0 + r)*(1.0 + r)*(-1.0 + s)*(1.0 + s)*t*(1.0 + t))/2., &
        -((-1.0 + r)*(1.0 + r)*(-1.0 + s)*(1.0 + s)*(-1.0 + t)*(1.0 + t)) /)
    Nr = (/ ((-1. + 2.*r)*(-1. + s)*s*(-1. + t)*t)/8.,              ((1. + 2.*r)*(-1. + s)*s*(-1. + t)*t)/8.,               &
        ((1. + 2.*r)*s*(1. + s)*(-1. + t)*t)/8.,                    ((-1. + 2.*r)*s*(1. + s)*(-1. + t)*t)/8.,               &
        ((-1. + 2.*r)*(-1. + s)*s*t*(1. + t))/8.,                   ((1. + 2.*r)*(-1. + s)*s*t*(1. + t))/8.,                &
        ((1. + 2.*r)*s*(1. + s)*t*(1. + t))/8.,                     ((-1. + 2.*r)*s*(1. + s)*t*(1. + t))/8.,                &
        -(r*(-1. + s)*s*(-1. + t)*t)/2.,                            ((1. - 2.*r)*(-1. + s)*(1. + s)*(-1. + t)*t)/4.,        &
        ((1. - 2.*r)*(-1. + s)*s*(-1. + t)*(1. + t))/4.,            -((1. + 2.*r)*(-1. + s**2)*(-1. + t)*t)/4.,             &
        -((1. + 2.*r)*(-1. + s)*s*(-1. + t**2))/4.,                 -(r*s*(1. + s)*(-1. + t)*t)/2.,                         &
        -((1. + 2.*r)*s*(1. + s)*(-1. + t**2))/4.,                  ((1. - 2.*r)*s*(1. + s)*(-1. + t)*(1. + t))/4.,         &
        -(r*(-1. + s)*s*t*(1. + t))/2.,                             ((1. - 2.*r)*(-1. + s)*(1. + s)*t*(1. + t))/4.,         &
        -((1. + 2.*r)*(-1. + s**2)*t*(1. + t))/4.,                  -(r*s*(1. + s)*t*(1. + t))/2.,                          &
        r*(-1. + s**2)*(-1. + t)*t,                                 r*(-1. + s)*s*(-1. + t**2),                             &
        ((-1. + 2.*r)*(-1. + s**2)*(-1. + t**2))/2.,                ((1. + 2.*r)*(-1. + s**2)*(-1. + t**2))/2.,             &
        r*s*(1. + s)*(-1. + t**2),                                  r*(-1. + s**2)*t*(1. + t),                              &
        -2.*r*(-1. + s**2)*(-1. + t**2) /)
    Ns = (/ ((-1. + r)*r*(-1. + 2.*s)*(-1. + t)*t)/8.,              (r*(1. + r)*(-1. + 2.*s)*(-1. + t)*t)/8.,               &
        (r*(1. + r)*(1. + 2.*s)*(-1. + t)*t)/8.,                    ((-1. + r)*r*(1. + 2.*s)*(-1. + t)*t)/8.,               &
        ((-1. + r)*r*(-1. + 2.*s)*t*(1. + t))/8.,                   (r*(1. + r)*(-1. + 2.*s)*t*(1. + t))/8.,                &
        (r*(1. + r)*(1. + 2.*s)*t*(1. + t))/8.,                     ((-1. + r)*r*(1. + 2.*s)*t*(1. + t))/8.,                &
        ((-1. + r)*(1. + r)*(1. - 2.*s)*(-1. + t)*t)/4.,            -((-1. + r)*r*s*(-1. + t)*t)/2.,                        &
        ((-1. + r)*r*(1. - 2.*s)*(-1. + t)*(1. + t))/4.,            -(r*(1. + r)*s*(-1. + t)*t)/2.,                         &
        (r*(1. + r)*(1. - 2.*s)*(-1. + t)*(1. + t))/4.,             -((-1. + r**2)*(1. + 2.*s)*(-1. + t)*t)/4.,             &
        -(r*(1. + r)*(1. + 2.*s)*(-1. + t**2))/4.,                  -((-1. + r)*r*(1. + 2.*s)*(-1. + t**2))/4.,             &
        ((-1. + r)*(1. + r)*(1. - 2.*s)*t*(1. + t))/4.,             -((-1. + r)*r*s*t*(1. + t))/2.,                         &
        -(r*(1. + r)*s*t*(1. + t))/2.,                              -((-1. + r**2)*(1. + 2.*s)*t*(1. + t))/4.,              &
        (-1. + r**2)*s*(-1. + t)*t,                                 ((-1. + r**2)*(-1. + 2.*s)*(-1. + t**2))/2.,            &
        (-1. + r)*r*s*(-1. + t**2),                                 r*(1. + r)*s*(-1. + t**2),                              &
        ((-1. + r**2)*(1. + 2.*s)*(-1. + t**2))/2.,                 (-1. + r**2)*s*t*(1. + t),                              &
        -2.*(-1. + r**2)*s*(-1. + t**2) /)
    Nt = (/((-1. + r)*r*(-1. + s)*s*(-1. + 2.*t))/8.,               (r*(1. + r)*(-1. + s)*s*(-1. + 2.*t))/8.,               &
        (r*(1. + r)*s*(1. + s)*(-1. + 2.*t))/8.,                    ((-1. + r)*r*s*(1. + s)*(-1. + 2.*t))/8.,               &
        ((-1. + r)*r*(-1. + s)*s*(1. + 2.*t))/8.,                   (r*(1. + r)*(-1. + s)*s*(1. + 2.*t))/8.,                &
        (r*(1. + r)*s*(1. + s)*(1. + 2.*t))/8.,                     ((-1. + r)*r*s*(1. + s)*(1. + 2.*t))/8.,                &
        ((-1. + r)*(1. + r)*(-1. + s)*s*(1. - 2.*t))/4.,            ((-1. + r)*r*(-1. + s)*(1. + s)*(1. - 2.*t))/4.,        &
        -((-1. + r)*r*(-1. + s)*s*t)/2.,                            (r*(1. + r)*(-1. + s)*(1. + s)*(1. - 2.*t))/4.,         &
        -(r*(1. + r)*(-1. + s)*s*t)/2.,                             ((-1. + r)*(1. + r)*s*(1. + s)*(1. - 2.*t))/4.,         &
        -(r*(1. + r)*s*(1. + s)*t)/2.,                              -((-1. + r)*r*s*(1. + s)*t)/2.,                         &
        -((-1. + r**2)*(-1. + s)*s*(1. + 2.*t))/4.,                 -((-1. + r)*r*(-1. + s**2)*(1. + 2.*t))/4.,             &
        -(r*(1. + r)*(-1. + s**2)*(1. + 2.*t))/4.,                  -((-1. + r**2)*s*(1. + s)*(1. + 2.*t))/4.,              &
        ((-1. + r**2)*(-1. + s**2)*(-1. + 2.*t))/2.,                (-1. + r**2)*(-1. + s)*s*t,                             &
        (-1. + r)*r*(-1. + s**2)*t,                                 r*(1. + r)*(-1. + s**2)*t,                              &
        (-1. + r**2)*s*(1. + s)*t,                                  ((-1. + r**2)*(-1. + s**2)*(1. + 2.*t))/2.,             &
        -2.*(-1. + r**2)*(-1. + s**2)*t /)
    
  end subroutine el12_ShapeFunc
    
    
subroutine ShapeFunc_Expand(NN, nen, N)
    ! expand [NN] -> [N]
    
    implicit none    
#include "Flash.h"    
    integer, intent(in) :: nen
    real, intent(in)  :: NN(nen)
    real, intent(out) :: N(NDIM, nen*NDIM)
    !integer :: i,j

    N = 0.
    ! New Numbering format:
    N(1,1:nen) = NN
    N(2,nen+1:2*nen) = NN
#if NDIM == 3	
    N(3,2*nen+1:nen*N_DIM) = NN
#endif
  end subroutine ShapeFunc_Expand
