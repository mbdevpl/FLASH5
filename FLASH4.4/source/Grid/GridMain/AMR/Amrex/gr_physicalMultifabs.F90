module gr_physicalMultifabs
    use amrex_amr_module, ONLY : amrex_multifab

    implicit none

    public :: unk
    public :: facevarx
    public :: facevary
    public :: facevarz

    type(amrex_multifab), allocatable, target :: unk(:)
    type(amrex_multifab), allocatable         :: facevarx(:)
    type(amrex_multifab), allocatable         :: facevary(:)
    type(amrex_multifab), allocatable         :: facevarz(:)
end module gr_physicalMultifabs

