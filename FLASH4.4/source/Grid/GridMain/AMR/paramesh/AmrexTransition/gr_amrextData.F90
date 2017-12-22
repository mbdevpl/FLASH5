module gr_amrextData
  use amrex_multifab_module, ONLY: amrex_multifab, amrex_multifab_destroy
  implicit none

  type(amrex_multifab),allocatable,TARGET :: gr_amrextUnkMFs(:)

contains

  subroutine gr_amrextDataInit(maxLev)
    integer,intent(IN) :: maxLev
    if (.not.allocated(gr_amrextUnkMFs)) then
       allocate(gr_amrextUnkMFs(0:maxLev-1))
    end if

  end subroutine gr_amrextDataInit

  subroutine gr_amrextDataFree()
    integer :: maxLev, level
    if (allocated(gr_amrextUnkMFs)) then
       maxLev = ubound(gr_amrextUnkMFs,1)
       do level = maxLev, 1, -1
          call amrex_multifab_destroy(gr_amrextUnkMFs(level-1))
       end do
       deallocate(gr_amrextUnkMFs)
    end if

  end subroutine gr_amrextDataFree
end module gr_amrextData
