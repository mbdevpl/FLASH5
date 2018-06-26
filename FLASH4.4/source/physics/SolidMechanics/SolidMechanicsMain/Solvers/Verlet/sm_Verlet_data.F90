
module sm_Verlet_data

! Structural Data per body:
type sm_Verlet_type

      !Timestep for structure:
      real :: dt

      ! Define integration variables (sub step info, etc)

end type sm_Verlet_type

! One to one correspondence with gr_sbBodyInfo for aero-grids.
type(sm_Verlet_type), save, dimension(:), pointer :: sm_Verlet_Info

end module sm_Verlet_data
