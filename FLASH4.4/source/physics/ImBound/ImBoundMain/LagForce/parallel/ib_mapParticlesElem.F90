

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

Subroutine ib_mapParticlesElem(eltype,vert_elem,nXi,nEta,ptelem,&
                              xi,yi,zi,ui,vi,wi,udi,vdi,wdi,     &
                              nxLi,nyLi,nzLi,                    &
                              xpos,ypos,zpos,xvel,yvel,zvel,     &
                              xacc,yacc,zacc,xnrm,ynrm,znrm,areai)


  implicit none

  integer, intent(IN)  :: eltype,vert_elem,nXi,nEta,ptelem
  real, intent(IN)     :: xi(vert_elem),yi(vert_elem),zi(vert_elem), &
                          ui(vert_elem),vi(vert_elem),wi(vert_elem), &
                          udi(vert_elem),vdi(vert_elem),wdi(vert_elem), &
                          nxLi(vert_elem),nyLi(vert_elem),nzLi(vert_elem)
  real, intent(OUT)   :: xpos(ptelem),ypos(ptelem),zpos(ptelem), &
                          xvel(ptelem),yvel(ptelem),zvel(ptelem), &
                          xacc(ptelem),yacc(ptelem),zacc(ptelem), &
                          xnrm(ptelem),ynrm(ptelem),znrm(ptelem), &
                          areai(ptelem)

  real :: area, area_sub, ei, ej, N1, N2, N3, N4

  integer :: i,ii,j,tmp

  real :: v12(MDIM), v13(MDIM), av(MDIM), indi, indj

  real :: nrm

  ! Select Case regarding Aerodynamic surface element type:
  select case(eltype)

#if NDIM == 2

  case(TWO_NODE_SEGMENT_VERT)

     ! The particles are defined in sub-segment vertices
     !   1    2   ...  ...  nXi   
     !   |----o----o----o----|     
     area = sqrt( (xi(2)-xi(1))**2. + (yi(2)-yi(1))**2.)

     do i = 1, nXi

        ei  = real(i-1)/real(nXi-1)

        ! Shape functions:
        N1 = 1. - ei
        N2 = ei
                 
        ! x,y,z positions of internal particles:
        xpos(i) = N1*xi(1) + N2*xi(2)
        ypos(i) = N1*yi(1) + N2*yi(2)
        zpos(i) = 0. 

        ! x,y,z velocities of internal particles:
        xvel(i) = N1*ui(1) + N2*ui(2) 
        yvel(i) = N1*vi(1) + N2*vi(2)
        zvel(i) = 0.

        ! x,y,z accelerations of internal particles:
        xacc(i) = N1*udi(1) + N2*udi(2) 
        yacc(i) = N1*vdi(1) + N2*vdi(2)
        zacc(i) = 0.        

        ! x,y,z normals of internal particles:
        xnrm(i) = N1*nxLi(1) + N2*nxLi(2) 
        ynrm(i) = N1*nyLi(1) + N2*nyLi(2)
        znrm(i) = 0.   


        ! Area associated with particle:
        areai(i) = area/real(nXi-1)
        if ((i .eq. 1).or.(i .eq. nXi)) areai(i) = areai(i)/2.

     enddo


  case(TWO_NODE_SEGMENT_CEN)

     ! The particles are defined in sub-segment centers
     !      1     2    ...   nXi   
     !   |--x--o--x--o--x--o--x--|     

     area = sqrt( (xi(2)-xi(1))**2. + (yi(2)-yi(1))**2.)

     do i = 1, nXi

        ei  = (0.5+real(i-1))/real(nXi)

        ! Shape functions:
        N1 = 1. - ei
        N2 = ei
                 
        ! x,y,z positions of internal particles:
        xpos(i) = N1*xi(1) + N2*xi(2)
        ypos(i) = N1*yi(1) + N2*yi(2)
        zpos(i) = 0. 

        ! x,y,z velocities of internal particles:
        xvel(i) = N1*ui(1) + N2*ui(2) 
        yvel(i) = N1*vi(1) + N2*vi(2)
        zvel(i) = 0.

        ! x,y,z accelerations of internal particles:
        xacc(i) = N1*udi(1) + N2*udi(2) 
        yacc(i) = N1*vdi(1) + N2*vdi(2)
        zacc(i) = 0.        

        ! x,y,z normals of internal particles:
        xnrm(i) = N1*nxLi(1) + N2*nxLi(2) 
        ynrm(i) = N1*nyLi(1) + N2*nyLi(2)
        znrm(i) = 0.   


        ! Area associated with particle:
        areai(i) = area/real(nXi)

     enddo


#elif NDIM == 3

  case(THREE_NODE_TRIANG_VERT)


     ! The particles are defined in sub-triang vertices
     !
     ! nEta  -
     !         \
     !       |  \
     !       o   o
     !             \
     !       |      \
     !       o   o   o
     !                 \
     !       |          \   
     !     2 o   o   o   o 
     !                     \
     !       |              \
     !     1 |---o---o---o---|  
     !       1   2     ...  nXi   



     ! Triangles Area
     v12(IAXIS) = xi(2)-xi(1);
     v12(JAXIS) = yi(2)-yi(1);
     v12(KAXIS) = zi(2)-zi(1);

     v13(IAXIS) = xi(3)-xi(1);
     v13(JAXIS) = yi(3)-yi(1);
     v13(KAXIS) = zi(3)-zi(1);

     call cross_pr(v12,v13,av)

     area = 0.5*sqrt(av(IAXIS)**2+av(JAXIS)**2+av(KAXIS)**2)

     area_sub = area/(real(nEta-1)**2.)

     ii = 0
     do j = 1, nEta
        indj = 0.
        if ((j .eq. 1) .or. (j .eq. nEta)) indj = 0.5

        tmp = 1 - j + nEta
        do i = 1, tmp

           ii = ii + 1
           
           ! Natural Coordinates
           ei  = real(i-1)/real(nEta-1)
           ej  = real(j-1)/real(nEta-1)

           ! Shape functions:
           N1 = 1. - ei - ej
           N2 = ei
           N3 = ej


           ! x,y,z positions of internal particles:
           xpos(ii) = N1*xi(1) + N2*xi(2) + N3*xi(3)
           ypos(ii) = N1*yi(1) + N2*yi(2) + N3*yi(3)
           zpos(ii) = N1*zi(1) + N2*zi(2) + N3*zi(3)

           ! x,y,z velocities of internal particles:
           xvel(ii) = N1*ui(1) + N2*ui(2) + N3*ui(3) 
           yvel(ii) = N1*vi(1) + N2*vi(2) + N3*vi(3)
           zvel(ii) = N1*wi(1) + N2*wi(2) + N3*wi(3)

           ! x,y,z accelerations of internal particles:
           xacc(ii) = N1*udi(1) + N2*udi(2) + N3*udi(3) 
           yacc(ii) = N1*vdi(1) + N2*vdi(2) + N3*vdi(3) 
           zacc(ii) = N1*wdi(1) + N2*wdi(2) + N3*wdi(3)

           ! x,y,z normals of internal particles:
           xnrm(ii) = N1*nxLi(1) + N2*nxLi(2) + N3*nxLi(3) 
           ynrm(ii) = N1*nyLi(1) + N2*nyLi(2) + N3*nyLi(3) 
           znrm(ii) = N1*nzLi(1) + N2*nzLi(2) + N3*nzLi(3)

           nrm = sqrt(xnrm(ii)**2 + ynrm(ii)**2 + znrm(ii)**2)

           xnrm(ii) = xnrm(ii) / nrm
           ynrm(ii) = ynrm(ii) / nrm
           znrm(ii) = znrm(ii) / nrm
           

           ! Area associated with particle:
           indi = 0.
           if ((i .eq. 1).or.(i .eq. tmp)) indi = 0.5
           areai(ii) = (2.-real(ceiling(indi+indj)))*area_sub
        end do
     end do

     areai(1)      = areai(1)/3.
     areai(nEta)   = areai(nEta)/3.
     areai(ptelem) = areai(ptelem)/3.



  case(THREE_NODE_TRIANG_CEN)

!!$     ! X2 - X1
!!$     Lseg = sqrt( (xi(2)-xi(1))**2. + (yi(2)-yi(1))**2. + (zi(2)-zi(1))**2.)
!!$     nXi  = ceiling(Lseg/dseg)
!!$     ! X3 - X1
!!$     Lseg = sqrt( (xi(3)-xi(1))**2. + (yi(3)-yi(1))**2. + (zi(3)-zi(1))**2.)
!!$     nEta = ceiling(Lseg/dseg)
!!$     nEta = max(nXi,nEta)
!!$     nXi  = nEta
!!$   
!!$     ! nXi is equal to nEta
!!$     ptelem = (nXi**2 + nXi)/2 


  case(FOUR_NODE_QUAD_VERT)

     call Driver_abortFlash("Four Node Quad Vertices, Not developed Yet")


  case(FOUR_NODE_QUAD_CEN)

     call Driver_abortFlash("Four Node Quad Cen, Not developed Yet")


#endif
  

  case default

     call Driver_abortFlash("Invalid Aerodynamic element type")

  endselect



  return

End Subroutine ib_mapParticlesElem

! -----------------------------------------------------------------------
! Function cross product:
! Computes cross produce of a(3) x b(3) 
!
! -----------------------------------------------------------------------
subroutine cross_pr(a,b,cr_pr)

  implicit none
  real a(3),b(3)
  real cr_pr(3)
  
  cr_pr(1) = a(2)*b(3) - a(3)*b(2)
  cr_pr(2) = a(3)*b(1) - a(1)*b(3)
  cr_pr(3) = a(1)*b(2) - a(2)*b(1)
      
  return
end subroutine cross_pr
