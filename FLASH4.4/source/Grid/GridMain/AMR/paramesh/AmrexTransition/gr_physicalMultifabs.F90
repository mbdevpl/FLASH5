!!****ih* source/Grid/GridMain/AMR/paramesh/AmrexTransition/gr_physicalMultifabs
!!
!! NAME
!!
!!  gr_physicalMultifabs
!!
!! SYNOPSIS
!!
!!  use gr_physicalMultifabs
!!
!! DESCRIPTION
!!
!!  A module that holds data structures for UNK and other global grid data
!!  structures, in their AMReX incarnation as  multifab arrays
!!
!! NOTES
!!
!!  UNK and other global grid data structures have traditionally been
!!  declared in a "physicaldata" module by FLASH, as multidimensional
!!  arrays of real data.
!!
!!  This implementation should be the same as the GridMain/AMR/Amrex
!!  implementation. Currently they differ in some TARGET attribute(s).
!!
!!***


module gr_physicalMultifabs
    use amrex_amr_module, ONLY : amrex_multifab

    implicit none

    public :: unk
    public :: facevarx
    public :: facevary
    public :: facevarz

    type(amrex_multifab), allocatable,TARGET :: unk(:)
    type(amrex_multifab), allocatable :: facevarx(:)
    type(amrex_multifab), allocatable :: facevary(:)
    type(amrex_multifab), allocatable :: facevarz(:)
end module gr_physicalMultifabs

