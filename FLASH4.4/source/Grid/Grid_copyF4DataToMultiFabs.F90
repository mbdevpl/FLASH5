subroutine Grid_copyF4DataToMultiFabs(gds, phi, nodetype, reverse)

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use tree,      ONLY : lrefine_max, lrefine

  use amrex_multifab_module

  type(amrex_multifab),intent(INOUT) :: phi(:)
  integer,intent(IN),OPTIONAL :: gds
  integer,intent(IN),OPTIONAL :: nodetype
  logical,intent(IN),OPTIONAL :: reverse

  ! Stub implementation, call is ignored!

end subroutine Grid_copyF4DataToMultiFabs
