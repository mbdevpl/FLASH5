module SolidMechanics_rbc_data 

implicit none
! some structural parameters 
real,parameter :: kc=2.4e-19;
real,parameter :: kBoltz=1.38e-23
real,parameter :: sm_rbc_gamma=4.5
real,parameter :: sm_rbc_Lambda=0.65
real,parameter :: sm_rbc_mass=1.;
real,parameter :: rho_w=993.68; ! kg/m^3
real,parameter :: sm_rbc_Y_M=392.453;
real,parameter :: sm_rbc_Y_P=18.9e-6;
real,parameter :: nu_mem_P=.022;


integer,parameter :: sm_rbc_Nf=27344
integer,parameter :: sm_rbc_m=2
integer,parameter :: sm_rbc_ka=4900,sm_rbc_kv=5000,sm_rbc_kd=100
integer,parameter :: sm_rbc_VR=94, sm_rbc_AR=135 !needed only in non-stress-free models
integer,parameter,dimension(3,3):: sm_rbc_Eye= RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))

real, save :: kbf, kbf_M, nu_mem_M;
real, save :: ks,fext, rho_Rbc, kBT
real, save :: r_M=7.82e-6/8.25, t_M, N_M, KBT_M, fext_M
logical, save:: Drag_flag =.false.
real, save :: drag_coeff = 0.;
integer, save :: m_exp=2;
real, parameter :: sm_rbc_T=296.  ! 
real, parameter :: mu_op=6.3e-6   ! Real units
real, save :: sm_rbc_gamma_t=90., sm_rbc_gamma_c=30., sm_rbc_sigma
real, save :: dpd_c1, dpd_c2, sm_rbc_mi
real, save :: mu_M, mu_p=3.2e-3 ;
real, save :: sm_rbc_nff=3., sm_rbc_nw=3.;
real, save :: sm_rbc_ae, sm_rbc_aw, sm_rbc_af=25.;
real, save :: sm_rbc_alpha_t=0.75, sm_rbc_viscCoeff=0.5;
integer, save :: sm_rbc_viscmodel=0
integer,save ::  rbcplotOutputInterval
integer,save ::  stretching_exp;
end module SolidMechanics_rbc_data
