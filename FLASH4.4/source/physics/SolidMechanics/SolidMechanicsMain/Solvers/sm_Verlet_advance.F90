! Stub

subroutine sm_Verlet_advance(ibd,restart)
      use Driver_interface, only : Driver_abortFlash
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart

      call Driver_AbortFlash("Verlet_advance not configured")

end subroutine sm_Verlet_advance

