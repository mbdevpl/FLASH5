!!
!! ib_countParticles:
!!
!!
!! NAME
!!  ib_countParticles
!!
!! SYNOPSIS
!!
!!  Counts the number of marker particles for one solid body.
!!
!!
!!*** 


#include "constants.h"
#include "Flash.h"

Subroutine ib_countParticles(ib,numPart)


  use gr_sbData, ONLY : gr_sbBodyInfo, NodesPerElem
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getDeltas
  use Grid_interface, ONLY : Grid_getMinCellSize
  use ib_interface, ONLY : ib_countParticlesElem

  implicit none

  integer, intent(IN)  :: ib
  integer, intent(OUT) :: numPart
  
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer :: count

  real :: del(MDIM)

  integer :: lb,blockID,aelem(NodesPerElem+CONSTANT_ONE)
  
  real :: Dmin

  integer :: vert_elem,nXi,nEta,ptelem

  real, dimension(NodesPerElem) :: xi,yi,zi

  integer :: i

  ! Get the Eulerian cell size
  call Grid_getMinCellSize(Dmin)


  ! We use the minimum eulerian cell size to define the number of points to be used on each 
  ! Aelem. Get Particle counts per element and particles per side Xi, Eta
  numPart = 0

  do i=1,gr_sbBodyInfo(ib)%NumAelem


     aelem(:) = gr_sbBodyInfo(ib)%AELEM(:,i)

     vert_elem = aelem(1)

     xi(1:vert_elem) = gr_sbBodyInfo(ib)%xb(aelem(2:vert_elem+1))
     yi(1:vert_elem) = gr_sbBodyInfo(ib)%yb(aelem(2:vert_elem+1))
#if NDIM == 3
     zi(1:vert_elem) = gr_sbBodyInfo(ib)%zb(aelem(2:vert_elem+1))
#endif 


     call ib_countParticlesElem(gr_sbBodyInfo(ib)%AELTYPE(i),vert_elem, &
                                xi(1:vert_elem),yi(1:vert_elem),zi(1:vert_elem),Dmin,nXi,nEta,ptelem)


     gr_sbBodyInfo(ib)%sbPtNumXi(i)  = nXi
     gr_sbBodyInfo(ib)%sbPtNumeta(i) = nEta
     gr_sbBodyInfo(ib)%sbPtNumElem(i)= ptelem

     numPart = numPart + ptelem

  enddo

  return
end Subroutine ib_countParticles
