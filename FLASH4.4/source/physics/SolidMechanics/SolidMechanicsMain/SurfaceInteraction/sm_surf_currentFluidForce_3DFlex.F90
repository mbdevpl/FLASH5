!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_currentFluidForce_3DFlex
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

subroutine sm_surf_currentFluidForce_3DFlex(ibd, ndofs, Hs_pres, Hs_visc)
  use SolidMechanics_data,         only: sm_bodyInfo, sm_nen, sm_structure
  use gr_sbData,                   ONLY: gr_sbBodyInfo
  use Driver_interface,            only: Driver_abortFlash
  use sm_element_interface,        only: sm_el10_FluidForce,sm_el02_FluidForce
  implicit none

#include "Flash_mpi.h"

  ! IO variables 
  integer, intent(in)  :: ibd, ndofs
  real,    intent(out) :: Hs_pres(ndofs), Hs_visc(ndofs)

  ! Local Variables
  integer :: p,e, nen_e, nee, a, i, pL, pG, neq
  real, allocatable, dimension(:) :: Hsp_pres, Hsp_visc
  real, dimension(NPART_PROPS) :: particle
  type(sm_structure), pointer :: body

  ! Meta-body gather buffer
  real, allocatable, dimension(:) :: sendbuf
  integer :: ierr

  ! assign body pointer
  body => sm_bodyInfo(ibd)

  ! clear the fluid storage
  Hs_pres = 0.
  Hs_visc = 0.

  ! Loop over particles
  do p = 1,gr_sbBodyInfo(ibd)%totalPart

     ! get patch info
     particle = gr_sbBodyInfo(ibd)% particles(1:NPART_PROPS,p)
     e        = int( particle(ELEM_PART_PROP) )
     nen_e    = sm_nen( body%ws_eltype(e) )
     nee      = body%max_dofs_per_node*nen_e  ! number of dofs per element, this is based on disp. dofs

     ! create storage for element force vector
     allocate( Hsp_pres(nee), Hsp_visc(nee) )

     select case( body%ws_eltype(e) )

#if NDIM == MDIM

     case( THREE_NODE_TRIANGLE  )           
        call sm_el02_FluidForce(body, e, particle, Hsp_pres, Hsp_visc )

     case( NINE_NODE_QUADRILATERAL )
        call sm_el10_FluidForce(body, e, particle, Hsp_pres, Hsp_visc )

#endif

     case default
        call Driver_abortFlash("WetSurface Element type not implemented.")

     end select

     ! Reduce Hsp into global force vector
     do a = 1,nen_e
        do i = 1,body%max_dofs_per_node

           ! local equation number 
           ! (groups dofs together as [u_1 u_2 ..., v_1 v_2 ..., w_1 w_2 ... ]
           !  this is how things are grouped everywhere in 3Dflexible)
           pL = a + nen_e*(i-1)

           ! Global Eqn Number
           pG = body%ID(i,body%ws_IEN(a,e))

           ! Store the fluid force for later
           Hs_pres(pG) = Hs_pres(pG) + Hsp_pres(pL)
           Hs_visc(pG) = Hs_visc(pG) + Hsp_visc(pL)

        end do
     end do

     ! Cleanup workspace
     deallocate( Hsp_pres, Hsp_visc )

  end do


  ! Here if this body is part of a metabody we will perform an all reduce sum for all 
  ! Pressure and viscous forces. 
  ! All processors managing bodies part of a metabody will run for their particular ibd
  ! through the routines sm_IntegX_advance, sm_assemble_ExtForce, sm_assembleFluidForces_rigid
  ! in this case. At this point they should have their contributions to body%Hs_pres, body%Hs_visc
  ! computed. An all reduce on the group is performed to gather the total corresponding
  ! fluid forces for the set:
  if (body%Metabody .gt. CONSTANT_ZERO) then
     allocate(sendbuf(body%ndofs))

     ! Pressure forces:
     sendbuf(:) = Hs_pres(:)
     call MPI_allReduce(sendbuf, Hs_pres, body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

     ! Viscous forces:
     sendbuf(:) = Hs_visc(:)
     call MPI_allReduce(sendbuf, Hs_visc, body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

     deallocate(sendbuf)
  endif

  return

end subroutine sm_surf_currentFluidForce_3DFlex
