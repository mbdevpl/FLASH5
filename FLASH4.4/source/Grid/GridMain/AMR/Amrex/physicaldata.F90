module physicaldata
    use amrex_amr_module, ONLY : amrex_multifab

    implicit none

    public :: unk

    type(amrex_multifab), allocatable :: unk(:)
end module physicaldata

