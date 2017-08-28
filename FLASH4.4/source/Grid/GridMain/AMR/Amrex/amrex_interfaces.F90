module amrex_interfaces
      interface
         subroutine gr_amrex_init()
            implicit none
         end subroutine gr_amrex_init
      end interface
      
      interface
         subroutine gr_amrex_finalize()
            implicit none
         end subroutine gr_amrex_finalize
      end interface
      
      interface
         subroutine gr_makeNewLevelFromScratch(lev, time, pba, pdm) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,     intent(IN), value :: lev
            real(wp),    intent(in), value :: time
            type(c_ptr), intent(in), value :: pba
            type(c_ptr), intent(in), value :: pdm
         end subroutine gr_makeNewLevelFromScratch 
      end interface
 
      interface
         subroutine gr_makeNewLevelFromCoarse(lev, time, pba, pdm) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,     intent(IN), value :: lev
            real(wp),    intent(in), value :: time
            type(c_ptr), intent(in), value :: pba
            type(c_ptr), intent(in), value :: pdm
         end subroutine gr_makeNewLevelFromCoarse 
      end interface
 
      interface
         subroutine gr_remakeLevel(lev, time, pba, pdm) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,     intent(IN), value :: lev
            real(wp),    intent(in), value :: time
            type(c_ptr), intent(in), value :: pba
            type(c_ptr), intent(in), value :: pdm
         end subroutine gr_remakeLevel 
      end interface 
 
      interface
         subroutine gr_clearLevel(lev) bind(c)
            implicit none
            integer, intent(in), value :: lev
         end subroutine gr_clearLevel
      end interface
 
      interface
         subroutine gr_markRefineDerefine(lev, tags, time, &
                                          tagval, clearval) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,           intent(IN), value :: lev
            type(c_ptr),       intent(in), value :: tags 
            real(wp),          intent(in), value :: time
            character(c_char), intent(in), value :: tagval
            character(c_char), intent(in), value :: clearval
         end subroutine gr_markRefineDerefine
      end interface
 
end module amrex_interfaces
