!!****if* source/Grid/GridMain/AMR/Amrex/gr_physicalMultifabs
!!
!! NAME
!!  gr_physicalMultifabs
!!
!! SYNOPSIS
!!
!!  use gr_physicalMultifabs
!!
!! DESCRIPTION 
!!  
!!  Define variables that store physical data as arrays of AMReX
!!  multifabs or as arrays of AMReX flux registers.
!!
!!  In all cases, the array index is the refinement level of the associated
!!  element.  Following AMReX, the level zero corresponds to the coarsest level
!!  of refinement.
!!
!!***

module gr_physicalMultifabs
    use amrex_amr_module,          ONLY : amrex_multifab
    use amrex_fluxregister_module, ONLY : amrex_fluxregister

    implicit none

    public :: unk
    public :: facevarx
    public :: facevary
    public :: facevarz

    type(amrex_multifab), allocatable, target :: unk(:)
    type(amrex_multifab), allocatable         :: facevarx(:)
    type(amrex_multifab), allocatable         :: facevary(:)
    type(amrex_multifab), allocatable         :: facevarz(:)

    type(amrex_fluxregister), allocatable :: flux_registers(:)
end module gr_physicalMultifabs

