

subroutine sm_io_rbcinit()
	
  use SolidMechanics_rbc_data, Only :  rbcplotOutputInterval, &
                                       stretching_exp,        &
                                       Fext

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  call RuntimeParameters_get("rbcplotOutputInterval", rbcplotOutputInterval)
  call RuntimeParameters_get("stretching_exp", stretching_exp)
  call RuntimeParameters_get("Fext", Fext) 
end subroutine sm_io_rbcinit
