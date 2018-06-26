! Stub

subroutine sm_Verlet_init(ibd,restart)
      use Driver_interface, only : Driver_abortFlash
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart

      call Driver_AbortFlash("Verlet_init not configured")

end subroutine sm_Verlet_init

