!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/2D/sm_el01_mapParticles
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
#include "constants.h"

subroutine sm_el01_mapParticles(body, e, ptelem,  &
                                xpos,ypos,        &
                                xvel,yvel,        &
                                xacc,yacc,        &
                                xnrm,ynrm,        &
                                areai, loc_num )
  use SolidMechanics_data, only: sm_structure
  use Driver_interface, only: Driver_abortFlash
  implicit none
  
  ! IO Variables
  type(sm_structure)   :: body     ! entire body structure
  integer, intent(in)  :: e        ! element number
  integer, intent(in)  :: ptelem
  real, dimension(ptelem) :: xpos, ypos, xvel, yvel, xacc, yacc, xnrm, ynrm, areai
  integer, dimension(ptelem) :: loc_num

  ! Internal variables
  integer, parameter :: nen_e = TWO_NODES
  integer :: a, i, idx1, idx2, nEta
  real, dimension(nen_e,NDIM) :: u, ud, udd
  real, dimension(NDIM) :: norm,v12
  real :: del,area,area_sub
  real :: ei,N1,N2

  !
  ! Get absolute position of each node in the element
  !
  do a = 1,nen_e
     idx1 = body%ws_IEN(a,e)
     u(a,1) = body%x(idx1)
     u(a,2) = body%y(idx1)
     do i = 1,NDIM
        idx2     = body%ID(i,idx1)
        u(a,i)   = u(a,i) + body%qn(idx2)
        ud(a,i)  = body%qdn(idx2)
        udd(a,i) = body%qddn(idx2)
     end do
  enddo

  ! nXi -> nEta
  nEta  = body%ws_nXi(e)
  
  !Get the normal:
  v12(IAXIS)=u(2,IAXIS)-u(1,IAXIS);
  v12(JAXIS)=u(2,JAXIS)-u(1,JAXIS);

  ! Segment Length
  area = sqrt( v12(IAXIS)**2. + v12(JAXIS)**2. )

  ! Normal outside: Counter clockwise node numbering in 2D 
  ! Therefore n = v12 x k 
  norm(IAXIS) =  v12(JAXIS)/area
  norm(JAXIS) = -v12(IAXIS)/area

  area_sub = area/real(ptelem);
  del=1./real(nEta);

  ! Map nEta particles:
  do i = 1, nEta

     ei  = (0.5+real(i-1))*del

     ! Shape functions:
     N1 = 1. - ei
     N2 = ei  

     ! x,y positions of internal particles:
     xpos(i) = N1*u(1,IAXIS) + N2*u(2,IAXIS);
     ypos(i) = N1*u(1,JAXIS) + N2*u(2,JAXIS);

     ! Vel
     xvel(i) = N1*ud(1,IAXIS) + N2*ud(2,IAXIS);
     yvel(i) = N1*ud(1,JAXIS) + N2*ud(2,JAXIS);

     ! Acc
     xacc(i) = N1*udd(1,IAXIS) + N2*udd(2,IAXIS);
     yacc(i) = N1*udd(1,JAXIS) + N2*udd(2,JAXIS);

     xnrm(i) = norm(IAXIS) ;
     ynrm(i) = norm(JAXIS) ;

     loc_num(i)=i;
     areai(i)  = area_sub;

  enddo

  return

end subroutine sm_el01_mapParticles
