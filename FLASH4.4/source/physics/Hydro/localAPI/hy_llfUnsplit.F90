!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_Unsplit/hy_llfUnsplit
!!
!! NAME
!!
!!  hy_llfUnsplit
!!
!! SYNOPSIS
!!
!!  call hy_hllUnsplit( integer (IN) :: blockCount,
!!                      integer (IN) :: blockList(blockCount),
!!                      real    (IN) :: dt,
!!                      real    (IN) :: dtOld  )
!!
!! DESCRIPTION
!!
!!  Performs Hydro update in a directionally unsplit fashion over a set
!!  of blocks.
!!  dt gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.
!! 
!!
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!  dtOld      -  old timestep (not used here)
!!
!! NOTES
!!
!!  This is a stb for a simple demo version. See documentation in the
!!  implementation file for more information.
!!
!! HISTORY
!!
!!  June  2013  - created KW
!!***


Subroutine hy_llfUnsplit ( blockCount, blockList, dt, dtOld )

  implicit none

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) ::  blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList
  real,    INTENT(IN) :: dt, dtOld
  !! -----------------------------------------------------


End Subroutine hy_llfUnsplit
