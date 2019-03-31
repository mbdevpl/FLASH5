!!****if* source/Grid/GridMain/AMR/Amrex/gr_initSpecific
!!
!! NAME
!!  gr_initSpecific
!!
!! SYNOPSIS
!!
!!  call gr_initSpecific()
!!
!!
!! DESCRIPTION
!!  Initialize some implementation-specific Grid data
!!
!!  This routine should initialize data in the gr_specificData module.
!!
!! ARGUMENTS
!!  none
!!
!!***


subroutine gr_initSpecific()
  use Grid_data, ONLY : gr_meshNumProcs
  use gr_specificData, ONLY : gr_nToLeft

  implicit none


  !Initialize grid arrays used by IO
  allocate(gr_nToLeft(0:gr_meshNumProcs-1))


end subroutine gr_initSpecific

