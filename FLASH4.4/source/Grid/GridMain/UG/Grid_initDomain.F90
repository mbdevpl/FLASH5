!!****if* source/Grid/GridMain/UG/Grid_initDomain
!!
!! NAME
!!
!!  Grid_initDomain
!!
!!
!! SYNOPSIS
!!
!!  call Grid_initDomain(logical(IN)  :: restart,
!!                       logical(INOUT) :: particlesInitialized)
!!
!!
!! DESCRIPTION
!! 
!!  Create the mesh, initialize all the mesh data structures
!!  and apply initial conditions
!!
!!
!! ARGUMENTS
!!
!!  restart : is true if the execution is starting from a checkpoint
!!            file, otherwise false.
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_initDomain( restart,particlesInitialized)

  use Grid_interface, ONLY : Grid_getLeafIterator, Grid_releaseLeafIterator, &
       Grid_fillGuardCells, Grid_getBlkPtr, Grid_releaseBlkPtr
  
  use Grid_data, ONLY : gr_gid, gr_eosModeInit, gr_eosMode, gr_eosModeNow, gr_globalMe,&
       gr_globalNumProcs
  use Eos_interface, ONLY : Eos_wrapped
  use Simulation_interface, ONLY : Simulation_initBlock, Simulation_initRestart
  use block_metadata, ONLY : block_metadata_t
  use leaf_iterator, ONLY : leaf_iterator_t
  implicit none

#include "Flash.h"
#include "constants.h"
  type(block_metadata_t) :: block
  type(leaf_iterator_t) :: itor
  
  logical, intent(in) :: restart
  logical, intent(inout) :: particlesInitialized

  integer :: blockID=1, ngid
  
  integer, dimension(LOW:HIGH, MDIM) :: blkLimits,blkLimitsGC
  real, pointer:: solnData(:,:,:,:)


  call gr_createDataTypes()

  if (NDIM == 1) then
        ngid = 5
     else if (NDIM == 2) then
        ngid = 9
     else if (NDIM == 3) then
        ngid = 15
  end if

  allocate(gr_gid(ngid,1)) !1 because in UG only 1 block per proc     

  if(.not.restart) then
     !  zero data in case we don't initialize all data in
     !  Simulation_initBlock... in particular the total vs. internal
     !  energies can cause problems in the eos call that follows
     call Grid_getLeafIterator(itor)
     call itor%blkMetaData(block)
     call Grid_getBlkPtr(block, solnData,CENTER)
     solnData = 0.0
     blkLimits(:,:) = block%limits
     call Simulation_initBlock(solnData,block)

     call Eos_wrapped(gr_eosModeInit, blkLimits, solnData)
     call Grid_releaseBlkPtr(block, solnData,CENTER)
     call Grid_releaseLeafIterator(itor)
  else
     ! Do no call the EOS here any more when restarting, since those calls
     ! may (depending on the Eos implementation, or whether Eos has
     ! been called explicitly before the checkpoint was written)
     ! introduce small data differences.
     ! As long as checkpoint data is written in a state where the
     ! solution data is thermodynamically consistent, there is no
     ! need for the call here. - KW

     ! Now give user code an opportunity to modify any values on restart.
     ! User implementation of Simulation_initRestart is responsible for
     ! calling Eos if that is needed to leave the solution data in a
     ! consistent state.

     call Simulation_initRestart()

  end if

  call Grid_fillGuardCells( CENTER_FACES, ALLDIR)

  gr_eosModeNow = gr_eosMode
  call gr_solversInit()


end subroutine Grid_initDomain
