!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/gr_sbCreateParticles
!!
!! NAME
!!  gr_sbCreateParticles
!!
!! SYNOPSIS
!!
!!  gr_sbCreateParticles()
!!  
!! DESCRIPTION 
!!  
!!  This routine is called from Grid_initDomain and Driver_evolveFlash for
!!  moving body. The master processor creates 
!!  particles for each triangle in solid body. Also gets the position coordinates in the
!!  particle structure. 
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

Subroutine gr_sbCreateParticles()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, aelem, totalPart,sumPart, &
                        gr_sbFirstCall, gr_sbStencil
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshMe, gr_meshComm

  use sm_surf_interface, only: sm_surf_countParticles_wsIEN, sm_surf_mapParticles, &
                               sm_surf_repoParticlesBC
  use SolidMechanics_data, only: sm_bodyInfo
  use sm_integinterface, only: sm_ga_isPredictor, sm_pc_isPredictor

  implicit none
  include "Flash_mpi.h"
  integer :: b,  flag_countParticles
  integer :: sumPart_pr, ierr

  call Timers_start("body_create_particles")

  do b = 1, gr_sbNumBodies

     ! IF fixed body and not the first call:
     if ((gr_sbBodyInfo(b)%sbIsFixed .eq. CONSTANT_ONE) .and. (gr_sbFirstCall .eq. CONSTANT_ZERO)) cycle

     !The solid body master processor creates the particles
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then

        ! Count the number of particles necessary for this body:

        flag_countParticles = SM_FALSE

        select case( sm_bodyInfo(b)%IntegMethod )
        case( SOLIDINTEG_GENALPHA )
           call sm_ga_isPredictor( b, flag_countParticles )
        case( SOLIDINTEG_PREDCORR )
           call sm_pc_isPredictor( b, flag_countParticles )
        case( SOLIDINTEG_MODVVERLET )
            flag_countParticles = SM_TRUE
        case default
           call Driver_abortFlash("Integrator type unknown.")
        end select 

        if( flag_countParticles == SM_TRUE  ) then

           call sm_surf_countParticles_wsIEN(b,totalPart)

           gr_sbBodyInfo(b) % totalPart = totalPart

           if (b .eq. 1) print*,gr_meshMe,"Total # of Particles: ",gr_sbBodyInfo(b)%totalPart
        
           ! Allocate the particles array:
           allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,gr_sbBodyInfo(b)%totalPart))

           if (gr_sbBodyInfo(b)%sbIsFixed .eq. CONSTANT_ONE) then
             allocate(gr_sbBodyInfo(b)%ielem(gr_sbStencil,MDIM,MDIM,gr_sbBodyInfo(b)%totalPart))
             allocate(gr_sbBodyInfo(b)%phile(gr_sbStencil,NDIM+1,MDIM,gr_sbBodyInfo(b)%totalPart))
             gr_sbBodyInfo(b)%ielem = NONEXISTENT
             gr_sbBodyInfo(b)%phile = 0.      
           endif

           gr_sbBodyInfo(b) % particles(:,:) = NONEXISTENT
           
           gr_sbBodyInfo(b) % particles(PEXT_PART_PROP,:) = 0. ! needed to make sure first corrector stop in GA

        end if

        ! Now do the Particle Mapping:
        call sm_surf_mapParticles(b)


        ! Test if Particles lie outside of domain, and reposition if there are periodic boundary conditions:
        call sm_surf_repoParticlesBC(b)

 
     else

        gr_sbBodyInfo(b) % totalPart = 0

     end if

  end do

  ! Get Sum of Particles for all bodies on Processor:
  sumPart_pr = 0
  do b = 1, gr_sbNumBodies
     sumPart_pr = sumPart_pr + gr_sbBodyInfo(b) % totalPart
  enddo
  call mpi_allreduce ( sumPart_pr, sumPart, 1, FLASH_INTEGER, &
                       MPI_MAX, gr_meshComm, ierr )

  call Timers_stop("body_create_particles")

  return

End Subroutine gr_sbCreateParticles
