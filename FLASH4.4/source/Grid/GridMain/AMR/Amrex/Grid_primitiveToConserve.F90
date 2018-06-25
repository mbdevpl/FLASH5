!!****if* source/Grid/GridMain/AMR/Amrex/Grid_primitiveToConserve
!!
!! NAME
!!  Grid_primitiveToConserve
!!
!!
!! SYNOPSIS
!!  Grid_primitiveToConserve(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: force)
!!
!!
!! DESCRIPTION
!!  Calls gr_primitiveToConserve
!!
!! ARGUMENTS
!!   blkList - integer list of blocks to be operated on
!!   count - number of blocks in the blkList
!!   force - whether to force conversion
!!
!!***

subroutine Grid_primitiveToConserve(blkList, count, force)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(IN) :: count
  integer, intent(IN) :: blkList(count)
  logical, intent(IN) :: force

  ! DEV: TODO This needs to be rethought or modernized to work with iterators
  call Driver_abortFlash("[Grid_primitiveToConserve] Not implemented for AMReX")
end subroutine Grid_primitiveToConserve

