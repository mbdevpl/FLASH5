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

    implicit none
    
    real(wp) :: t_old, t_new

    call cpu_time(t_old)
    write(*,*) "Holy moly! It runs!"
    call cpu_time (t_new)
    write (*,*) ' Cpu time = ', (t_new - t_old)
end subroutine Driver_evolveFlash

