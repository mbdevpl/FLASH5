! DEV: TODO This should be renamed.  It looks like it comes from the AMReX
!           as is.
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
         subroutine gr_initNewLevelCallback(lev, time, pba, pdm) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,     intent(IN), value :: lev
            real(wp),    intent(in), value :: time
            type(c_ptr), intent(in), value :: pba
            type(c_ptr), intent(in), value :: pdm
         end subroutine gr_initNewLevelCallback
      end interface
 
      interface
         subroutine gr_makeFineLevelFromCoarseCallback(lev, time, pba, pdm) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,     intent(IN), value :: lev
            real(wp),    intent(in), value :: time
            type(c_ptr), intent(in), value :: pba
            type(c_ptr), intent(in), value :: pdm
         end subroutine gr_makeFineLevelFromCoarseCallback
      end interface
 
      interface
         subroutine gr_remakeLevelCallback(lev, time, pba, pdm) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,     intent(IN), value :: lev
            real(wp),    intent(in), value :: time
            type(c_ptr), intent(in), value :: pba
            type(c_ptr), intent(in), value :: pdm
         end subroutine gr_remakeLevelCallback
      end interface 
 
      interface
         subroutine gr_clearLevelCallback(lev) bind(c)
            implicit none
            integer, intent(in), value :: lev
         end subroutine gr_clearLevelCallback
      end interface
 
      interface
         subroutine gr_markRefineDerefineCallback(lev, tags, time, &
                                          tagval, clearval) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            integer,           intent(IN), value :: lev
            type(c_ptr),       intent(in), value :: tags 
            real(wp),          intent(in), value :: time
            character(c_char), intent(in), value :: tagval
            character(c_char), intent(in), value :: clearval
         end subroutine gr_markRefineDerefineCallback
      end interface

      interface
         subroutine gr_fillPhysicalBC(pmf, scomp, ncomp, time, pgeom) bind(c)
            use iso_c_binding
            use amrex_fort_module, ONLY : wp => amrex_real
            implicit none
            type(c_ptr),    value :: pmf
            type(c_ptr),    value :: pgeom
            integer(c_int), value :: scomp
            integer(c_int), value :: ncomp
            real(wp),       value :: time
         end subroutine gr_fillPhysicalBC
      end interface
 
      interface
         subroutine gr_getFinestLevel(level)
            implicit none
            integer, intent(IN) :: level
         end subroutine gr_getFinestLevel
      end interface

      interface
         subroutine gr_writeData(stepno, t_new)
            implicit none
            integer, intent(IN) :: stepno
            real,    intent(IN) :: t_new
         end subroutine gr_writeData
      end interface

end module amrex_interfaces

