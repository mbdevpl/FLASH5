!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/sm_el02_FluidForce
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

subroutine sm_el02_FluidForce(body, e, particle, Hsp_pres, Hsp_visc )
  use SolidMechanics_data, only: sm_structure

  implicit none

  ! constants
  integer, parameter :: nen_e = THREE_NODES
  integer, parameter :: nee   = NDIM*nen_e

  ! IO Variables
  type(sm_structure),  intent(in) :: body
  integer, intent(in) :: e  ! element number on the wet surface
  real, dimension(NPART_PROPS), intent(in)  :: particle
  real, dimension(nee), intent(out) :: Hsp_pres, Hsp_visc

  ! Internal Variables
  integer :: a, i, j,k, idx1, idx2, i1, p, nxi, neta
  real, dimension(nen_e,NDIM) :: u
  real :: N1,N2,N3,del,delp,ei,ej
  real, dimension(NDIM) ::  N, nhat, f_pres, f_visc
  real, external :: ddot, dnrm2
  integer :: pL,tmp,ii
  !
  ! Get absoluate position of each node in the element
  !
  do a = 1,nen_e
     idx1 = body%ws_IEN(a,e)
     
     do i = 1,3
        idx2     = body%ID(i,idx1)
        u(a,i)   = body%qn(idx2)
     end do
  enddo

  ! Get patch parameters
  nxi  = body%ws_nXi(e)
  neta = body%ws_nEta(e)
  
  ! local patch number
  p = int( particle(PLOC_PART_PROP) ) 


  ! Get the shape functions for this element.
  del=1./real(nEta);
  ii = 0
  do j=1,nEta
     tmp = 1 - j + nEta;
     do i = 1,tmp
        
        ii = ii + 1;
        if (ii.eq.p) then
           !write(*,*) 'That is the poiny',ii,p
           delp=(1.-1.5*del)/real(nEta-1);
           
           ei=0.5*del+ real(i-1) *delp;
           ej=0.5*del+ real(j-1) *delp;
           ! Shape functions:
           N1 = 1. - ei - ej;
           N2 = ei;
           N3 = ej;
        end if
     end do
  end do
  N=(/N1,N2,N3/)
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
   !write(*,*)' f_pres(a)',f_pres
   !write(*,*)' f_visc(a)', f_visc
   
   do i1 = 1,nen_e
      
      ! Pressure forces
      Hsp_pres(NDIM*i1-2) = N(i1) * f_pres(IAXIS) * particle(AREA_PART_PROP);
      Hsp_pres(NDIM*i1-1) = N(i1) * f_pres(JAXIS) * particle(AREA_PART_PROP);
      Hsp_pres(NDIM*i1  ) = N(i1) * f_pres(KAXIS) * particle(AREA_PART_PROP);
      ! Viscous forces
      Hsp_visc(NDIM*i1-2) = N(i1) * f_visc(IAXIS) * particle(AREA_PART_PROP);
      Hsp_visc(NDIM*i1-1) = N(i1) * f_visc(JAXIS) * particle(AREA_PART_PROP);
      Hsp_visc(NDIM*i1  ) = N(i1) * f_visc(KAXIS) * particle(AREA_PART_PROP);
      
   end do
   
 end subroutine sm_el02_FluidForce
!!$        do a = 1,NDIM
!!$        pL = a + nen_e*(i1-1)
!!$        write(*,*) 'pL',pL,  f_pres(a),f_visc(a)
!!$        Hsp_pres(pL) = Hsp_pres(pL) + N(i1) * f_pres(a) * particle(AREA_PART_PROP);
!!$        Hsp_visc(pL) = Hsp_visc(pL) + N(i1) * f_visc(a) * particle(AREA_PART_PROP);
!!$        
!!$     end do
