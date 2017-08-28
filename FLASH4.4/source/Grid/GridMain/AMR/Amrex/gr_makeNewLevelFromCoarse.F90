subroutine gr_makeNewLevelFromCoarse(lev, time, ba, dm) bind(c)
   use iso_c_binding
   use amrex_fort_module, ONLY : wp => amrex_real
   
   implicit none
   
   integer,     intent(IN), value :: lev
   real(wp),    intent(in), value :: time
   type(c_ptr), intent(in), value :: ba
   type(c_ptr), intent(in), value :: dm

   write(*,*) "[gr_makeNewLevelFromCoarse]"
end subroutine gr_makeNewLevelFromCoarse

