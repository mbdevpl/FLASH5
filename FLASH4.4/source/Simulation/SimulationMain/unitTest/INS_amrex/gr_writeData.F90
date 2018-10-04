#include "Flash.h"

subroutine gr_writeData(stepno, t_new, argBaseName)
    use amrex_fort_module,     ONLY : wp => amrex_real
    use amrex_amrcore_module,  ONLY : amrex_get_numlevels, &
                                      amrex_geom, &
                                      amrex_ref_ratio
    use amrex_string_module,   ONLY : amrex_string_build, &
                                      amrex_string
    use amrex_plotfile_module, ONLY : amrex_write_plotfile

    use gr_physicalMultifabs,  ONLY : unk

    implicit none

    integer,  intent(IN) :: stepno
    real(wp), intent(IN) :: t_new
    character(len=*), intent(IN), optional :: argBaseName

    character(17), parameter :: PLOT_FILE = "ins_amrex_plot_"

    integer              :: nlevs
    character(len=127)   :: filename
    character(len=16)    :: current_step
    character(4)         :: current_var

    type(amrex_string), allocatable :: varname(:)
    integer,            allocatable :: stepno_arr(:)

    integer :: i 

    if      (stepno .lt. 1000000) then
       write(current_step,fmt='(i5.5)') stepno
    else if (stepno .lt. 10000000) then
       write(current_step,fmt='(i6.6)') stepno
    else if (stepno .lt. 100000000) then
       write(current_step,fmt='(i7.7)') stepno
    else if (stepno .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') stepno
    else
       write(current_step,fmt='(i15.15)') stepno
    end if
    filename = trim(plot_file) // current_step

    nlevs = amrex_get_numlevels()

    allocate(varname(NUNK_VARS))

    do i = 1, SIZE(varname)
        write(current_var,'(I4.4)') i
        call amrex_string_build(varname(i), "var"//TRIM(current_var))
    end do

    allocate(stepno_arr(0:nlevs-1))
    stepno_arr = stepno

    call amrex_write_plotfile(filename, nlevs, unk, varname, amrex_geom, &
                              t_new, stepno_arr, amrex_ref_ratio)

    deallocate(varname)
    deallocate(stepno_arr)
end subroutine gr_writeData

