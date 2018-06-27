!
!
!
!!***

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

Subroutine ib_countParticlesElem(eltype,vert_elem,xi,yi,zi,Dmin,nXi,nEta,ptelem)

  use Driver_interface, ONLY : Driver_abortFlash

  use ImBound_Data, ONLY : distFactor

  implicit none

  integer, intent(IN)  :: eltype,vert_elem
  real, intent(IN)     :: xi(vert_elem),yi(vert_elem),zi(vert_elem),Dmin
  integer, intent(OUT) :: nXi,nEta,ptelem

  real :: dseg, Lseg
  
  dseg = distFactor*Dmin

  ! Select Case regarding Aerodynamic surface element type:
  select case(eltype)

#if NDIM == 2

  case(TWO_NODE_SEGMENT_VERT)

     ! The particles are defined in sub-segment vertices
     !   1    2   ...  ...  nXi   
     !   |----o----o----o----|     
     Lseg = sqrt( (xi(2)-xi(1))**2. + (yi(2)-yi(1))**2.)

     nXi  = ceiling(Lseg/dseg) + 1 ! Points along Aero segment
     nEta = nXi                   ! Not Used

     ptelem = nXi

  case(TWO_NODE_SEGMENT_CEN)

     ! The particles are defined in sub-segment centers
     !      1     2    ...   nXi   
     !   |--x--o--x--o--x--o--x--|     
     Lseg = sqrt( (xi(2)-xi(1))**2. + (yi(2)-yi(1))**2.)

     nXi  = ceiling(Lseg/dseg)     ! Points along Aero segment
     nEta = nXi                   ! Not Used

     ptelem = nXi



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

     ! X2 - X1
     Lseg = sqrt( (xi(2)-xi(1))**2. + (yi(2)-yi(1))**2. + (zi(2)-zi(1))**2.)
     nXi  = ceiling(Lseg/dseg) + 1
     ! X3 - X1
     Lseg = sqrt( (xi(3)-xi(1))**2. + (yi(3)-yi(1))**2. + (zi(3)-zi(1))**2.)
     nEta = ceiling(Lseg/dseg) + 1
     nEta = max(nXi,nEta)
     nXi  = nEta
   
     ! nXi is equal to nEta
     ptelem = (nXi**2 + nXi)/2


  case(THREE_NODE_TRIANG_CEN)

     ! X2 - X1
     Lseg = sqrt( (xi(2)-xi(1))**2. + (yi(2)-yi(1))**2. + (zi(2)-zi(1))**2.)
     nXi  = ceiling(Lseg/dseg)
     ! X3 - X1
     Lseg = sqrt( (xi(3)-xi(1))**2. + (yi(3)-yi(1))**2. + (zi(3)-zi(1))**2.)
     nEta = ceiling(Lseg/dseg)
     nEta = max(nXi,nEta)
     nXi  = nEta
   
     ! nXi is equal to nEta
     ptelem = (nXi**2 + nXi)/2 


  case(FOUR_NODE_QUAD_VERT)

     call Driver_abortFlash("Four Node Quad Vertices, Not developed Yet")


  case(FOUR_NODE_QUAD_CEN)

     call Driver_abortFlash("Four Node Quad Cen, Not developed Yet")


#endif
  



  case default

     call Driver_abortFlash("Invalid Aerodynamic element type")

  endselect

  return
End Subroutine ib_countParticlesElem
