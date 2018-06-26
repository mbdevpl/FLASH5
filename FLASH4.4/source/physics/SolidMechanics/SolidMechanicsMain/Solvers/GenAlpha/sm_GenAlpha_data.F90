
module sm_GenAlpha_data

! Structural Data per body:
type sm_GenAlpha_type

   ! Timestep for structure:
   real :: dt

   ! Spectral Radius
   real :: rhoinf

   ! solver parameters
   real :: gamma, beta, alpha_m, alpha_f

   ! Define integration variables
   real    :: pcerr,pcepsilon
   integer :: pcflag
   integer :: pcconvflag 

   ! Predictor/Corrector GA tolerance
   real    :: pred_epsilon, corr_epsilon

   ! relax the predictor
   real    :: pred_pos_relax

   ! override incase full fsi isn't needed 
   integer :: predictor_only

   ! override incase body is fixed during some transient calculations
   integer :: zeroQn_only

end type sm_GenAlpha_type

! One to one correspondence with gr_sbBodyInfo for aero-grids.
type(sm_GenAlpha_type), save, dimension(:), pointer :: sm_GenAlpha_Info

end module sm_GenAlpha_data
