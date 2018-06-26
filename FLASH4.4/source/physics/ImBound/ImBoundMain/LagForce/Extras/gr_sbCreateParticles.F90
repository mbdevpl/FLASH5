!!****if* source/physics/ImBound/ImBoundMain/LagForce/Extras/gr_sbCreateParticles
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

Subroutine gr_sbCreateParticles()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, aelem, totalPart,sumPart
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use ib_interface, ONLY : ib_countParticles, ib_mapParticles

  implicit none
  include "Flash_mpi.h"
  real :: xpos, ypos, zpos
  integer, save :: numParticles
  integer :: b, i, j, p
  logical, parameter :: localLogFile = .true.

  integer :: sumPart_pr,ierr

  integer :: tmp, numTriangles, ie 
  real :: ei, ej, N1, N2, N3
  integer, parameter :: nSideParticles = 2

  real :: xvel,yvel,zvel,area,areai
  integer :: ael_1,ael_2

  integer, save :: pt_maxPerProc
  
  call Timers_start("body_create_particles")

  call RuntimeParameters_get('pt_maxPerProc', pt_maxPerProc)

  do b = 1, gr_sbNumBodies

     !The solid body master processor creates the particles
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then

!==================================================================================================
!============================    KPD Change    ====================================================
!==================================================================================================
!gr_sbBodyInfo(1)%xo = 0.0
!gr_sbBodyInfo(1)%yo = 0.0
!gr_sbBodyInfo(1)%zo = 0.0
!print*,"Xo,Yo,Zo Defined AGAIN in gr_sbCreateParticles.F90 for core#",gr_sbBodyInfo(b) % bodyMaster
print*,"Xo,Yo,Zo NOT Defined again in gr_sbCreateParticles.F90",gr_sbBodyInfo(1)%xo,gr_sbBodyInfo(1)%yo
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        ! Count the number of particles necessary for this body:
        call ib_countParticles(b,totalPart)
        gr_sbBodyInfo(b) % totalPart = totalPart

        print*,gr_meshMe,"KPD # of Particles:",gr_sbBodyInfo(b)%totalPart,"ALLOC",pt_maxPerProc

        ! Allocate the particles array:
!==================================================================================================
!============================    KPD Change    ====================================================
!==================================================================================================
!       allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,gr_sbBodyInfo(b)%totalPart))
        allocate(gr_sbBodyInfo(b) % particles(1:NPART_PROPS,1:pt_maxPerProc))
        !                                                   ^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        gr_sbBodyInfo(b) % particles(:,:) = NONEXISTENT

        ! Now do the Particle Mapping:
        call ib_mapParticles(b)
 
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
