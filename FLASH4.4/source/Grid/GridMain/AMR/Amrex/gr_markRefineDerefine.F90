subroutine gr_markRefineDerefine(lev, tags, time, tagval, clearval) bind(c)
   use iso_c_binding
   use amrex_fort_module, ONLY : wp => amrex_real
   
   implicit none
   
   integer,           intent(IN), value :: lev
   type(c_ptr),       intent(in), value :: tags 
   real(wp),          intent(in), value :: time
   character(c_char), intent(in), value :: tagval
   character(c_char), intent(in), value :: clearval

   write(*,*) "[gr_markRefineDerefine]"
end subroutine gr_markRefineDerefine

