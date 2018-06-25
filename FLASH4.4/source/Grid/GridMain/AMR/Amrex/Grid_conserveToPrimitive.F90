!!****if* source/Grid/GridMain/AMR/Amrex/Grid_conserveToPrimitive
!!
!! NAME
!!  Grid_conserveToPrimitive
!!
!! SYNOPSIS
!!  Grid_conserveToPrimitive(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: allCells,
!!                           logical(in) :: force)
!!
!! DESCRIPTION
!!  Calls gr_conserveToPrimitive
!!
!! ARGUMENTS
!!   blkList - integer list of blocks to be operated on
!!   count - number of blocks in the blkList
!!   allCells - act on all cells, including guardcells, if .TRUE.,
!!              otherwise only modify interior cells.
!!   force - whether to force conversion
!!
!!***

subroutine Grid_conserveToPrimitive(blkList, count, allCells, force)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
  integer, intent(IN) :: count
  integer, intent(IN) :: blkList(count)
  logical, intent(IN) :: allCells
  logical, intent(IN) :: force

  ! DEV: TODO This needs to be rethought or modernized to work with iterators
  call Driver_abortFlash("[Grid_conserveToPrimitive] Not implemented for AMReX")
end subroutine Grid_conserveToPrimitive

