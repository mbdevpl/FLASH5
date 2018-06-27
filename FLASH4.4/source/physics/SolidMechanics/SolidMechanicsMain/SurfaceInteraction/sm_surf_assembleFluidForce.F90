!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_assembleFluidForce
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

subroutine sm_surf_assembleFluidForce(ibd)
  use SolidMechanics_data,         only: sm_bodyInfo, sm_nen, sm_structure
  use gr_sbData,                   ONLY: gr_sbBodyInfo
  use Driver_interface,            only: Driver_abortFlash
  use sm_element_interface,        only: sm_el10_FluidForce,sm_el02_FluidForce
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use sm_integdata,                only: sm_integ_subiter
  implicit none

#include "Flash_mpi.h"

  ! IO variables 
  integer, intent(in) :: ibd

  ! Local Variables
  integer :: p,e, nen_e, nee, a, i, pL, pG, neq
  real, allocatable, dimension(:) :: Hsp_pres, Hsp_visc
  real, dimension(NPART_PROPS) :: particle
  type(sm_structure), pointer :: body
  LOGICAL :: restart
  LOGICAL, save :: first_call = .true.

  ! Under-relaxation parameters
  real :: relax_param, smash
  real, parameter :: smash_50p = 10.
  real, parameter :: smash_stretch = 5.  
  real, parameter :: smash_scale = 1.0/(1.0 - tanh( (1. - smash_50p)/smash_stretch  ) )

  ! Meta-body gather buffer
  real, allocatable, dimension(:) :: sendbuf
  integer :: ierr

  ! assign body pointer
  body => sm_bodyInfo(ibd)
  neq = body%neq

  call RuntimeParameters_get("restart",restart)

  ! Don't do anything if there is no particles
  !if( gr_sbBodyInfo(ibd)%totalPart == 0 .and. restart ) return

  ! Get the relaxtion parameter
  call RuntimeParameters_get("relax_fluidforce_param",relax_param)

  ! Smash the relaxation parameter as the subiter count increases  
  ! smash = 1/2*(1-tanh( (n - 50p)/stretch )
  smash = smash_scale*(1.0 - tanh( ( real(sm_integ_subiter+1) - smash_50p)/smash_stretch  ) )
  relax_param = smash*relax_param

  ! Save old values of the fluid forces
  body%Hsi_pres = body%Hs_pres
  body%Hsi_visc = body%Hs_visc

  if( first_call .and. restart ) then
     first_call = .false.
     body%Hs(1:neq) =  body%Hs(1:neq) + body%Hs_visc(1:neq) + body%Hs_pres(1:neq)

  else

     ! clear the fluid storage
     body%Hs_pres = 0.
     body%Hs_visc = 0.
     write(*,*)'gr_sbBodyInfo(ibd)%totalPart',gr_sbBodyInfo(ibd)%totalPart
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

#else /* 2D */
        case( TWO_NODE_LINE )
           call Driver_abortFlash("This is not 2D safe...so far.")
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
              body%Hs_pres(pG) = body%Hs_pres(pG) + Hsp_pres(pL)
              body%Hs_visc(pG) = body%Hs_visc(pG) + Hsp_visc(pL)

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
        allocate(sendbuf(Body%ndofs))

        ! Pressure forces:
        sendbuf(:) = body%Hs_pres(:)
        call MPI_allReduce(sendbuf, body%Hs_pres, Body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

        ! Viscous forces:
        sendbuf(:) = body%Hs_visc(:)
        call MPI_allReduce(sendbuf, body%Hs_visc, Body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

        deallocate(sendbuf)
     endif


     ! Update the External force vector
     body%Hs(1:neq) =  body%Hs(1:neq) + relax_param*(      body%Hs_pres(1:neq)  + body%Hs_visc(1:neq)  ) &
          + (1.-relax_param)*( body%Hsi_pres(1:neq) + body%Hsi_visc(1:neq) )

  end if

#ifdef DEBUG_SOLID
#if NDIM == MDIM
  ! Write FX, FY, FZ
  write(*,'(A25,3g18.10)') 'Pressure Forces', &
       sum(body%Hs_pres(body%ID(IAXIS,:))),   &
       sum(body%Hs_pres(body%ID(JAXIS,:))),   &
       sum(body%Hs_pres(body%ID(KAXIS,:)))

  write(*,'(A25,3g18.10)') 'Viscous Forces', &
       sum(body%Hs_visc(body%ID(IAXIS,:))),  &
       sum(body%Hs_visc(body%ID(JAXIS,:))),  &
       sum(body%Hs_visc(body%ID(KAXIS,:)))
#else
  ! Write FX, FY
  write(*,'(A25,3g18.10)') 'Pressure Forces', &
       sum(body%Hs_pres(body%ID(IAXIS,:))),   &
       sum(body%Hs_pres(body%ID(JAXIS,:)))

  write(*,'(A25,3g18.10)') 'Viscous Forces', & 
       sum(body%Hs_visc(body%ID(IAXIS,:))),  &
       sum(body%Hs_visc(body%ID(JAXIS,:)))
#endif

  ! write out current Forcing relaxation parameter
  write(*,'(A25,g18.10)') 'Relaxation parameter =', relax_param

  !write(*,'(A25,g18.10,g18.10,A2)') 'Hs [',maxval(body%Hs),minval(body%Hs),']'

#endif

  return

end subroutine sm_surf_assembleFluidForce
