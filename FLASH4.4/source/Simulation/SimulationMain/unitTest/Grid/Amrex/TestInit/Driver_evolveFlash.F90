!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestIterator/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash ()
!!
!! DESCRIPTION
!!
!!  Simple stripped down version for testing single units.
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : wp => amrex_real
    use amrex_amr_module,      ONLY : amrex_geom, &
                                      amrex_problo, &
                                      amrex_probhi
    use amrex_amrcore_module,  ONLY : amrex_get_finest_level
    use amrex_parallel_module, ONLY : amrex_parallel_myproc, &
                                      amrex_parallel_nprocs

    implicit none
    
    integer  :: n_procs = -1 
    integer  :: rank = -1

    real(wp) :: t_old, t_new

    call cpu_time(t_old)
    write(*,*) "Holy moly! It runs!"

    n_procs = amrex_parallel_nprocs()
    rank = amrex_parallel_myproc()
    if (rank == 0) then
        write(*,*) 
        write(*,*) "Geometry information"
        write(*,*) "Domain Lower-Left  = ", amrex_problo
        write(*,*) "Domain Upper-Right = ", amrex_probhi
    end if

    call cpu_time (t_new)
    write (*,*) ' Cpu time = ', (t_new - t_old)
end subroutine Driver_evolveFlash

