module gr_amrextData
  use amrex_multifab_module, ONLY: amrex_multifab
  implicit none

  type(amrex_multifab),allocatable,TARGET :: gr_amrextUnkMFs(:)

contains

  subroutine gr_amrextDataInit(maxLev)
    integer,intent(IN) :: maxLev
    if (.not.allocated(gr_amrextUnkMFs)) then
       allocate(gr_amrextUnkMFs(maxLev))
    end if

  end subroutine gr_amrextDataInit
end module gr_amrextData
