!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/sm_el02_mapParticles
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

subroutine sm_el02_mapParticles(body, e, ptelem,  &
                                xpos,ypos,zpos,   &
                                xvel,yvel,zvel,   &
                                xacc,yacc,zacc,   &
                                xnrm,ynrm,znrm,   &
                                areai, loc_num )
  use SolidMechanics_data, only: sm_structure
  use Driver_interface, only: Driver_abortFlash
  use sm_misc_interface,only : sm_crossprod
  implicit none
  
  ! IO Variables
  type(sm_structure)   :: body     ! entire body structure
  integer, intent(in)  :: e        ! element number
  integer, intent(in)  :: ptelem
  real, dimension(ptelem) :: xpos, ypos, zpos, xvel, yvel, zvel, xacc, yacc, zacc, xnrm, ynrm, znrm, areai
  integer, dimension(ptelem) :: loc_num

  ! Internal variables
  integer, parameter :: nen_e = THREE_NODES
  integer :: a, i, ii, j, idx1, idx2, nEta, tmp
  real, dimension(nen_e,NDIM) :: u, ud, udd
  real, dimension(NDIM) :: norm,v12,v13
  real :: del,delp,area,area_sub
  real :: ei,ej,N1,N2,N3
  !
  ! Get absoluate position of each node in the element
  !
  do a = 1,nen_e
     idx1 = body%ws_IEN(a,e)
     u(a,1) = body%x(idx1)
     u(a,2) = body%y(idx1)
     u(a,3) = body%z(idx1)
     
     if (Body%BodyType.eq.BODYTYPE_RBC) u(a,1:NDIM)=0.;

     do i = 1,NDIM
        idx2     = body%ID(i,idx1)
        u(a,i)   = u(a,i) + body%qn(idx2)
        ud(a,i)  = body%qdn(idx2)
        udd(a,i) = body%qddn(idx2)
     end do
  enddo

  ! nXi = nEta
  nEta  = body%ws_nXi(e)
  
  !Get the normals
  v12(1)=u(2,IAXIS)-u(1,IAXIS);
  v12(2)=u(2,JAXIS)-u(1,JAXIS);
  v12(3)=u(2,KAXIS)-u(1,KAXIS);
  v13(1)=u(3,IAXIS)-u(1,IAXIS);
  v13(2)=u(3,JAXIS)-u(1,JAXIS);
  v13(3)=u(3,KAXIS)-u(1,KAXIS);

  
  norm=sm_crossprod(v12,v13);
  
  area=0.5*sqrt(sum(norm*norm));
  norm=norm/(2.*area);

  area_sub = area/real(ptelem);
  del=1./real(nEta);
  ii = 0
  do j=1,nEta
     tmp = 1 - j + nEta;
     do i = 1,tmp
        
        ii = ii + 1;
        delp=(1.-1.5*del)/real(nEta-1);
        
        ei=0.5*del+ real(i-1) *delp;
        ej=0.5*del+ real(j-1) *delp;
        ! Shape functions:
        N1 = 1. - ei - ej;
        N2 = ei;
        N3 = ej;
        
        
        ! x,y,z positions of internal particles:
        xpos(ii) = N1*u(1,IAXIS) + N2*u(2,IAXIS) + N3*u(3,IAXIS);
        ypos(ii) = N1*u(1,JAXIS) + N2*u(2,JAXIS) + N3*u(3,JAXIS) ;
        zpos(ii) = N1*u(1,KAXIS) + N2*u(2,KAXIS) + N3*u(3,KAXIS) ;

        ! Vel
        xvel(ii) = N1*ud(1,IAXIS) + N2*ud(2,IAXIS) + N3*ud(3,IAXIS);
        yvel(ii) = N1*ud(1,JAXIS) + N2*ud(2,JAXIS) + N3*ud(3,JAXIS) ;
        zvel(ii) = N1*ud(1,KAXIS) + N2*ud(2,KAXIS) + N3*ud(3,KAXIS) ;

        ! Acc
        
        xacc(ii) = N1*udd(1,IAXIS) + N2*udd(2,IAXIS) + N3*udd(3,IAXIS);
        yacc(ii) = N1*udd(1,JAXIS) + N2*udd(2,JAXIS) + N3*udd(3,JAXIS) ;
        zacc(ii) = N1*udd(1,KAXIS) + N2*udd(2,KAXIS) + N3*udd(3,KAXIS) ;


        xnrm(ii) = norm(IAXIS) ;
        ynrm(ii) = norm(JAXIS) ;
        znrm(ii) = norm(KAXIS) ; 
        loc_num(ii)=ii ;
        areai(ii) = area_sub;
     end do
  end do

  return

end subroutine sm_el02_mapParticles
