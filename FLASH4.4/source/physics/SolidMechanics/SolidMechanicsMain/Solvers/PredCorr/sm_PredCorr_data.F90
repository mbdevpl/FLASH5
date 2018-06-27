
module sm_PredCorr_data

! Structural Data per body:
type sm_PredCorr_type

      ! Timestep for structure:
      real :: dt

      ! Define integration variables
      ! Predictor-Corrector methods:
      integer :: pcmethod, pcflag, pciter
      real    :: pcerr,pcepsilon
      integer :: pcconvflag 
      real, allocatable, dimension(:) :: vardt

end type sm_PredCorr_type

! One to one correspondence with gr_sbBodyInfo for aero-grids.
type(sm_PredCorr_type), save, dimension(:), pointer :: sm_PredCorr_Info

end module sm_PredCorr_data
