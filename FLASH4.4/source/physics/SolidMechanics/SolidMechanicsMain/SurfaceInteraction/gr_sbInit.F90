!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/MovingBodies/gr_sbInit
!!
!! NAME
!!  gr_sbInit
!!
!! SYNOPSIS
!!
!!  call gr_sbInit()
!!
!! DESCRIPTION
!!
!!  * Called from Grid_initDomain
!!  * Read input file containing the coordinates of the vertices and the triangle elements of one body (ref body).
!!  * Get the min and max values of physical domain, and calculate the positions where the bodies need to be placed (in matrix form) 
!!  * Create triangles representing each body
!!  * Calculate centroids of each triangle
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbInit()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
                        gr_sbDebug,gr_sbParticleCount,gr_sbIJKtoProc,gr_sbStencil
  use SolidMechanics_data, only : sm_BodyInfo
  use ImBound_data, only : ib_stencil
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

#ifdef FLASH_GRID_UG  
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs, gr_axisMe, gr_axisNumProcs
#else
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
#endif

  implicit none
  include "Flash_mpi.h"
  
  integer :: b

  integer, allocatable, dimension(:,:,:) :: gr_sbIJKtoProcLoc 
  integer :: procs,ierr

  call RuntimeParameters_get("sm_NumBodies", gr_sbNumBodies)
  allocate(gr_sbBodyInfo(gr_sbNumBodies))
  allocate(gr_sbParticleCount(gr_sbNumBodies))
  ! Particle count to zero:
  gr_sbParticleCount(:) = 0

  ! First Set the number of vertices and triangles:
  do b = 1, gr_sbNumBodies
     gr_sbBodyInfo(b)%Mype = gr_meshMe
     gr_sbBodyInfo(b)%bodyMaster = sm_BodyInfo(b)%bodyMaster 

     ! Set value of IB interpolation stencil for allocation of ielem, phile
     ! Fields of gr_sbBodyInfo for fixed bodies. 
     gr_sbStencil = ib_stencil

     ! Set gr_sbBodyInfo(b)sbIsFixed
     gr_sbBodyInfo(b)%sbIsFixed = CONSTANT_ZERO  !!! NOW SET TO MOVING BODY, to
                                                 !!! be assigned from sm_BodyInfo(b)%IsFixed
     !gr_sbBodyInfo(b)%sbIsFixed = CONSTANT_ONE

     ! Dummy allocation to deallocate in gr_sbFinalize - Imbound.F90
     gr_sbBodyInfo(b)%totalpart = 0
     allocate(gr_sbBodyInfo(b)%particles(0,0))

     !if (gr_sbBodyInfo(b)sbIsFixed .eq. CONSTANT_ONE)) then
     !   allocate(bodyInfo%ielem(0,0,0,0))
     !   allocate(bodyInfo%phile(0,0,0,0))
     !endif


  enddo

#ifdef FLASH_GRID_UG  

  procs = gr_axisNumProcs(IAXIS)*gr_axisNumProcs(JAXIS)*gr_axisNumProcs(KAXIS)
  if (gr_meshNumProcs .ne. procs) call Driver_AbortFlash("gr_sbInit : for UG, procs not equal to gr_meshNumProcs")

  ! Construct array ijk to proc, mapping from ijk processor location to
  ! list of processors on the mesh. 
  allocate(gr_sbIJKtoProcLoc(gr_axisNumProcs(IAXIS),gr_axisNumProcs(JAXIS),gr_axisNumProcs(KAXIS)))
  allocate(gr_sbIJKtoProc   (gr_axisNumProcs(IAXIS),gr_axisNumProcs(JAXIS),gr_axisNumProcs(KAXIS)))

  gr_sbIJKtoProc    = 0
  gr_sbIJKtoProcLoc = 0
  gr_sbIJKtoProcLoc(gr_axisMe(IAXIS)+1,gr_axisMe(JAXIS)+1,gr_axisMe(KAXIS)+1) = gr_meshMe

  ! Allreduce so all procs know the other procs locations:
  call mpi_allreduce ( gr_sbIJKtoProcLoc, gr_sbIJKtoProc, procs, &
                       FLASH_INTEGER, MPI_SUM, gr_meshComm, ierr )

  deallocate(gr_sbIJKtoProcLoc)

#endif
  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
