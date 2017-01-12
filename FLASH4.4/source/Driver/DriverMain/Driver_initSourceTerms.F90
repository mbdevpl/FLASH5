!!****if* source/Driver/DriverMain/Driver_initSourceTerms
!!
!! NAME
!!   
!!  Driver_initSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_initSourceTerms(logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes all source terms Units by 
!!   calling their respective initialization routines
!!  
!! ARGUMENTS
!!   myPE - current processor number
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_initSourceTerms( restart)

  use Polytrope_interface, ONLY:  Polytrope_init
  use Burn_interface, ONLY:  Burn_init
  use Heat_interface, ONLY : Heat_init
  use Diffuse_interface, ONLY : Diffuse_init
  use Ionize_interface, ONLY : Ionize_init
  use Flame_interface, ONLY : Flame_init
  use RadTrans_interface, ONLY : RadTrans_init
  use Deleptonize_interface, ONLY : Deleptonize_init

  implicit none
  
  logical, intent(in) :: restart

  call Polytrope_init()
  call Diffuse_init()
  call Heat_init()
  call Ionize_init()
  call Burn_init()
  call Flame_init()
  call RadTrans_init()
  call Deleptonize_init()

end subroutine Driver_initSourceTerms
