!!****if* source/Grid/GridMain/paramesh/Grid_sendOutputData
!!
!! NAME
!!  Grid_sendOutputData
!!
!! SYNOPSIS
!!
!!  call Grid_sendOutputData()
!!
!! DESCRIPTION
!!
!!  This routine allows the Grid unit to checkpoint any scalar data
!!  stored in the Grid_data Fortran modules and is called by the
!!  routine IO_updateScalars before checkpointing.  To send data to
!!  the IO unit this routine calls IO_setScalar.  In addition this
!!  routine may prepare any other data the IO unit needs for
!!  checkpointing which the Grid unit owns.
!!
!!  For example, the Grid unit owns the variable which defines the
!!  grid geometry gr_geometry. This value needs to be checkpointed so
!!  that visualization tools can determine if the simulation ran with
!!  cartesian, spherical, cylindrical etc. coordinates.  To send
!!  scalar data such as gr_geometry to the IO unit to be checkpointed,
!!  the routine calls IO_setScalar("geometry", gr_geometry)
!!
!!  Other data which is sent to the IO unit is globalNumBlocks, globalOffset
!!  and globalNumBlocks
!!
!!  In implementation for PARAMESH, This routine also prepares the gid data structure which is also
!!  needed for visualization purposes. gid((nfaces + nchildren + 1 parent)
!!
!!
!!  ARGUMENTS
!!   none
!!
!!  SEE ALSO
!!
!!   IO_setScalar, IO_updateScalars
!!
!!***

subroutine Grid_sendOutputData()

#include "constants.h"
#include "Flash.h"

  use Grid_data, ONLY :gr_str_geometry, gr_globalNumBlocks, &
        gr_meshComm, gr_meshMe, gr_meshNumProcs
  use Grid_interface, ONLY : Grid_getLocalNumBlks

  use gr_specificData, ONLY : gr_nToLeft,gr_globalOffset
  use gr_specificData, ONLY : gr_ioLocalNumBlocks

  use IO_interface, ONLY : IO_setScalar

  implicit none

 include "Flash_mpi.h"

  integer :: localNumBlocks, i, ierr

  !Put NXB, NYB and NZB into a saved variable to prevent problems with
  !an xlf "feature."
  integer,parameter :: local_nxb = NXB
  integer,parameter :: local_nyb = NYB
  integer,parameter :: local_nzb = NZB
  integer,parameter :: dimensionality = NDIM


  !set the scalars for the grid unit
  call IO_setScalar("nxb", local_nxb)
  call IO_setScalar("nyb", local_nyb)
  call IO_setScalar("nzb", local_nzb)

  call IO_setScalar("geometry", gr_str_geometry)
  call IO_setScalar("dimensionality", dimensionality)


  call gr_updateData()          ! This (re)computes gr_ioLocalNumBlocks and gr_ioBlkBLAH arrays.

  !! Get the local number of blocks from everybody
  !! to find total number of blocks
  localNumBlocks = gr_ioLocalNumBlocks

  call MPI_Allgather(localNumBlocks, 1, MPI_INTEGER, gr_nToLeft, &
       1, MPI_INTEGER, gr_meshComm, ierr)

  gr_globalNumBlocks = 0

  do i = 0,gr_meshNumProcs-1
     gr_globalNumBlocks = gr_globalNumBlocks + gr_nToLeft(i)
  end do


  call IO_setScalar("globalNumBlocks", gr_globalNumBlocks)

  ! compute the number of processors to the left of a processor
  do i = gr_meshNumProcs-1,1,-1
     gr_nToLeft(i) = gr_nToLeft(i-1)
  end do


  gr_nToLeft(0) = 0
  do i = 2,gr_meshNumProcs-1
     gr_nToLeft(i) = gr_nToLeft(i) + gr_nToLeft(i-1)
  end do

  gr_globalOffset = gr_nToLeft(gr_meshMe)



  !need to make sure this happens before checkpointing
  call gr_ptFillBlkParticleInfo()


end subroutine Grid_sendOutputData
