subroutine gr_writeData(stepno, t_new)
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

    character(17), parameter :: PLOT_FILE = "sedov_amrex_plot_"

    integer              :: nlevs
    character(len=127)   :: filename
    character(len=16)    :: current_step
    type(amrex_string)   :: varname(4)
    integer, allocatable :: stepno_arr(:)

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

    call amrex_string_build(varname(1), "dens")
    call amrex_string_build(varname(2), "ener")
    call amrex_string_build(varname(3), "pres")
    call amrex_string_build(varname(4), "temp")

    allocate(stepno_arr(0:nlevs-1))
    stepno_arr = stepno

    call amrex_write_plotfile(filename, nlevs, unk, varname, amrex_geom, &
                              t_new, stepno_arr, amrex_ref_ratio)

end subroutine gr_writeData

