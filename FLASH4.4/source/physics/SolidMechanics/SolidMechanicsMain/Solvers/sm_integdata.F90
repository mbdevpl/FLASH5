
module sm_integdata
  ! Data for global and synchronization of different time integrators (which
  ! are a function of each body).


  ! Flag that defines if all bodies have been iterated to convergence on
  ! a given timestep:
  integer, save :: sm_convflag_all

  ! Max error among all bodies on a given sub-iteration
  real, save :: sm_errmax_all

  ! Number of Subiteration
  integer, save :: sm_integ_subiter
  integer, save :: sm_integ_subiter_old

  ! Error threshold to assume solution divergence and stop the run
  real, parameter :: sm_err_diverge = 1.e12 

end module sm_integdata
