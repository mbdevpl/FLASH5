!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/gr_sbFinalize
!!        source/physics/SolidMechanics/SolidMechanicsMain/SurfInteraction/gr_sbFinalize
!!
!! NAME
!!
!!  gr_sbFinalize
!!
!! SYNOPSIS
!!
!!  gr_sbFinalize()
!!
!! DESCRIPTION
!!
!!  Called from Grid_finalize. Deallocates the data structure that holds all information
!!  about the solid bodies.
!!
!! ARGUMENTS
!!
!!
!!***
#include "SolidMechanics.h"
#include "Flash.h"

Subroutine gr_sbFinalize()
  use gr_sbData,           ONLY: gr_sbBodyInfo, gr_sbNumBodies, gr_sbParticleCount, gr_sbFirstCall
  use SolidMechanics_data, only: sm_bodyInfo
  use sm_integinterface,   only: sm_ga_isPredictor, sm_pc_isPredictor

  implicit none
  include "Flash_mpi.h"
  integer :: b, ierr, flag_dealloc_particles

  do b = 1, gr_sbNumBodies

     ! IF fixed body and not the first call:
     if ((gr_sbBodyInfo(b)%sbIsFixed .eq. CONSTANT_ONE) .and. (gr_sbFirstCall .eq. CONSTANT_ZERO)) cycle
 

     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then

        flag_dealloc_particles = SM_FALSE

        select case( sm_bodyInfo(b)%IntegMethod )
        case( SOLIDINTEG_GENALPHA )
           call sm_ga_isPredictor( b, flag_dealloc_particles )

        case( SOLIDINTEG_PREDCORR )
           call sm_pc_isPredictor( b, flag_dealloc_particles )

        case( SOLIDINTEG_MODVVERLET )
            flag_dealloc_particles = SM_TRUE

        case default
           call Driver_abortFlash("gr_sbFinalize: Integrator type unknown.")
        end select       

        if( flag_dealloc_particles == SM_TRUE  ) then

           deallocate(gr_sbBodyInfo(b) % particles)     

        end if

     else
        if(gr_sbParticleCount(b) > 0) then
           deallocate(gr_sbBodyInfo(b) % particles)
        end if
     end if

  end do
  !if (allocated(gr_sbParticleCount)) deallocate(gr_sbParticleCount)
  !deallocate(gr_sbBodyInfo)

End Subroutine gr_sbFinalize
