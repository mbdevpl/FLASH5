!!****if* source/Simulation/SimulationMain/Sedov/sim_scaleProfile
!!
!! NAME
!!
!!  sim_scaleProfile
!!
!! SYNOPSIS
!!
!!  call sim_scaleProfile(real(in) :: tcurr)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   tcurr : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine sim_scaleProfile(tcurr)
  use Simulation_data, ONLY : sim_nProfile, &
       sim_useProfileFromFile,              &
       sim_profileIsScaled,sim_profileScaledTime, &
       sim_rhoAmbient, sim_expEnergy,       &
       sim_profileInitial,                  &
       sim_rProf, sim_vProf, sim_rhoProf, sim_pProf
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none

  real,intent(in)    :: tcurr

  real     ::  S, St            !scaling for analytical solution from file
  integer :: i

  ! Scale data from the solution file to the current time

  if (sim_useProfileFromFile) then
     if (.NOT. sim_profileIsScaled .OR. (sim_profileScaledTime .NE. tcurr)) then
        call Timers_start("scaleProfile")
        S  = ( sim_rhoAmbient / (sim_expEnergy*tcurr*tcurr) ) ** 0.2
        St = S * tcurr

        do i=1,sim_nProfile+1
           sim_rProf  (i) = sim_profileInitial(1,i) / S
           sim_vProf  (i) = sim_profileInitial(2,i) / St
           sim_rhoProf(i) = sim_profileInitial(3,i) * sim_rhoAmbient
           sim_pProf  (i) = sim_profileInitial(4,i) * sim_rhoAmbient / (St*St)
!!$        print*,i,sim_rProf  (i),sim_vProf  (i),sim_rhoProf(i),sim_pProf  (i)
!!$        if (.NOT. written) write(53,*)sim_rProf  (i),sim_vProf  (i),sim_rhoProf(i),sim_pProf  (i)
        end do
        call Timers_stop("scaleProfile")
     end if
  end if

  sim_profileIsScaled   = .TRUE.
  sim_profileScaledTime = tcurr

end subroutine sim_scaleProfile
