!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/SolidMechanics_computeDt
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!  
!!
!!***
#include "constants.h"
#include "Flash.h"  
#include "SolidMechanics.h"
  
subroutine SolidMechanics_computeDt(dt_solid)

  use SolidMechanics_data, only: sm_BodyInfo, sm_structure, sm_NumBodies, sm_MeshMe, sm_meshComm
  use sm_integinterface,   only: sm_EstDT_3DFlexible, sm_EstDT_rigid, sm_EstDT_rbc
  use Driver_interface,    only: Driver_getNStep
  implicit none

#include "Flash_mpi.h"

  ! IO Variables
  real, INTENT(INOUT) :: dt_solid

  ! Internal Variables
  integer :: ibd, nstep, ierr
  real    :: dt, dt_min

  ! Init all the variables
  dt = 1.e12
  dt_min = 1.e12
  
  ! get NSTEP
  call Driver_getNStep( nstep )
!!$  write(*,*) 'Calling the Solid compute DT'
!!$  write(*,*) 'nstep',nstep,'mod( nstep, 100 )',mod( nstep, 100 )
!!$  stop
!!$  if( mod( nstep, 100 ) == 0 ) then

     ! Loop over bodies:
     do ibd=1,sm_NumBodies 

        if ( sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster ) then

           if (sm_BodyInfo(ibd)%neq .gt. 0) then

           select case( sm_bodyInfo(ibd)%BodyType )

           case( BODYTYPE_RIGID ) 
              call sm_EstDT_rigid(ibd, dt)

#if NDIM == MDIM

           case( BODYTYPE_3DFLEXIBLE )
              call sm_EstDT_3DFlexible(ibd, dt)

           case( BODYTYPE_RBC )
              call sm_EstDT_rbc(ibd, dt)

#else /* 2D */

#endif
              
           case default
              call Driver_abortFlash("BodyType not yet implemented. <SolidMechanics_computeDt>")

           end select    

           dt_min = min( dt, dt_min )

           endif

        end if ! if proc has a body

     end do ! loop over bodies

     ! Find Minimum across procs
     call mpi_allreduce ( dt_min, dt_solid, 1, FLASH_REAL, &
                          MPI_MIN, sm_meshComm, ierr )
     
!!  end if ! if 'check' is being performed

  return

end subroutine SolidMechanics_computeDt
