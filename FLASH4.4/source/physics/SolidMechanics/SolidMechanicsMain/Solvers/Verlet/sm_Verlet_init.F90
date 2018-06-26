!

subroutine sm_Verlet_init(ibd,restart)
      use SolidMechanics_data, only: sm_NumBodies
      use sm_Verlet_data, only: sm_Verlet_type, sm_Verlet_info
      use Driver_interface, ONLY : Driver_getDt
      use Driver_data, Only : dr_dtInit
      implicit none
      ! Argument list
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
      ! local variables
      type(sm_Verlet_type), pointer :: integ
      real :: dt
      ! Check if sm_Verlet_info has been "allocated"
      if( .not. associated( sm_Verlet_info ) ) then
         allocate( sm_Verlet_info( sm_NumBodies ) )
      end if

      ! set things specific for your body here...
      integ => sm_Verlet_info(ibd)
      !call Driver_getDt(dt)
      !write(*,*) ' dt ',dt
      !write(*,'(A,f12.8,A)',advance='no') 'Initial Solid timestep= ',integ%dt,' ...    '
      !write(*,*) ' dt ',dt
      !dt=min(integ%dt,dr_dtInit)
      
      !integ%dt =dt; 
      !dt=integ%dt;
       dt=0.001;
       integ%dt=dt;
      write(*,'(A,f12.8)') 'Final Solid timestep= ',integ%dt

end subroutine sm_Verlet_init

