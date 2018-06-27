subroutine sm_rbc_Calc_parameters(ibd )
#define PRINTOUT_PARMS  
#include "SolidMechanics.h"
#include "constants.h"
#include "Flash.h"

  use SolidMechanics_data  
  use SolidMechanics_rbc_data 
  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_data, Only : dr_dt
  implicit none
  ! Argument list
  integer, intent(in) :: ibd
  
  
  ! Local variables 
  integer :: Nnodes, Nele
  integer ::i
  real :: kb_conc,Ka_conc,kv_conc,LL,Lm,Kp2
  real :: xo,Lo, Lmax, Dist_avg
  real :: conv_factor
  real :: x1,x2,y1,y2,z1,z2, vec(3)
  integer :: sp_modeltype
  integer :: node1, node2
  real, allocatable, dimension(:):: normals, centers, vol_centers,new_comm
  
  real :: comm(NDIM)
  type(sm_structure), pointer :: body
  
  ! pointer to sm_BodyInfo
  body => sm_BodyInfo(ibd)
  
  
  Nnodes= body%nnp;
  Nele  = body%nele;
  allocate(normals(Nele*NDIM))
  allocate(centers(Nele*NDIM))
  allocate(vol_centers(Nele*NDIM))
  allocate(new_comm(Nele*NDIM))
  
  body%tho =  acos((sqrt(3.0)*(NNodes-2.)-5.*PI)  &
       /(sqrt(3.0)*(NNodes-2.)-3.*PI)) !radians
 
  body%cos_tho=cos(Body%tho)
  body%sin_tho=sin(Body%tho)
  
  body%r_mult=2.05
  conv_factor=sqrt(real(sm_rbc_Nf-2.0)/real(NEle-2.0))
  
  kbf=(2./sqrt(3.0))*kc
  
  !
  !     Scaling basis
  !  
  ! Energy scaling base.
  kBT_M=sm_rbc_Y_P*(r_M * r_M)/sm_rbc_Y_M  
  ! Force scaling base.
  N_M=kBT_M/r_M;  
  ! Bending coefficient in DPD units 
  kbf_M=kbf/kBT_M;    
  ! Kb T in DPD units
  kBT=kBoltz*sm_rbc_T/kBT_M; 
  ! Density of rbc
  rho_RBC = 1.15 * rho_w ;
  ! mass of one RBC bead
  sm_rbc_mi = rho_RBC * sm_rbc_VR * 1.e-18 / NNodes; ! Kg
  ! Membrane viscosity of RBCs
  nu_mem_M=sqrt(3.)*(sm_rbc_gamma_t+(sm_rbc_gamma_c/4.));
  
  ! Time scaling bases
  t_M=((r_M*sm_rbc_Y_M*nu_mem_P)/( sm_rbc_Y_P*nu_mem_M))**sm_rbc_alpha_t;
  
  !  Dimensionless shear modulus 
  mu_M = mu_p * (r_M * r_M)/(N_M * t_M);
  if (Drag_flag) then
     drag_coeff = 6 * PI * mu_M ;
  else    
     drag_coeff=0;
  end if
  
  dpd_c1=sqrt(4.0*sm_rbc_gamma_t*KBT);
  
  dpd_c2=sqrt(2.0*KBT*(3.0*sm_rbc_gamma_c-sm_rbc_gamma_t));
  
  sm_rbc_sigma=sqrt(2.*sm_rbc_gamma*KBT);
  
  ! Amplitude of the conservative forces between water particles
  !af=25.;                    
  sm_rbc_ae=0.39*((sm_rbc_nff*kBT_M)+(0.1*sm_rbc_af*sm_rbc_nf*sm_rbc_nf))/  &
       ((0.0303*sm_rbc_nw**sm_rbc_nw)+(0.5617**sm_rbc_nw)-0.8536);
  sm_rbc_aw=(sm_rbc_ae*sm_rbc_ae)/sm_rbc_af;

  !     Squaring the cutoffs
!!$  RCUTSQ = RCUT * RCUT ;
!!$  LCUTSQ = LCUT * LCUT ;
  
  
  
#ifdef PRINTOUT_PARMS
  write(*,*) ' body%tho',body%tho
  write(*,*) ' body%cos_tho', body%cos_tho
  write(*,*) ' body%sin_tho', body%sin_tho
  write(*,*)
  write(*,*), 'Bending coefficient kbf= ',kbf
  write(*,*)' '
  write(*,*) 'Energy scaling basis KBT_M= ',KBT_M
  write(*,*)' '
  write(*,*) 'N_M ', N_M
  write(*,*)' '
  write(*,*)'The bending coeff kb_M ',kbf_M
  write(*,*)' '
  write(*,*)'KBT ',KBT
  write(*,*)' '
  write(*,*)'The mass scaling basis: ',sm_rbc_mi
  write(*,*)' '
  write(*,*)'The dissp coeff gamma_t: ',sm_rbc_gamma_t
  write(*,*)' '      
  write(*,*)'The dissp coeff gamma_C: ',sm_rbc_gamma_C
  write(*,*)' '
  write(*,*)'The scaling exponent alpha_t: ',sm_rbc_alpha_t
  write(*,*)' '
  write(*,*)'Model Membrane viscosity: ',nu_mem_M
  write(*,*)' '
  write(*,*)'Physical membrane viscosity: ',nu_mem_P
  write(*,*)' '
  write(*,*)'Time scaling basis: ',t_M
  write(*,*)' '
  write(*,*)'Actual time step: ',dr_dt*t_M,' s'
  write(*,*)''
  write(*,*)'The flow viscosity is: ',mu_p*1e3,'mPa.s'
  write(*,*)''
  write(*,*)'Fluid dynamic viscosity mu_M: ',mu_M
  write(*,*)' '
  write(*,*)'Drag coefficient: ',Drag_coeff
  write(*,*)' '
  write(*,*)' '
  write(*,*)'Random forces 1st coeff: ', dpd_c1
  write(*,*)' '
  write(*,*)' '
  write(*,*)'Random forces 2nd coeff: ', dpd_c2
  
#endif 
  
  
  !
  !     Spring model types
  !
  sp_modeltype = body%MatType(Nele+1);
  
  select case (sp_modeltype)
     
  case (MATERIAL_SF_FENE_POW , MATERIAL_SF_WLC_POW , MATERIAL_SF_WLC_FENE ,  MATERIAL_SF_WLC_CQ )

     ! Some calculations for the allocation;
     do i=1,Body%nLinks
        
        node1=body%Links(1,i)
        node2=body%Links(2,i)
        
        x1 = Body%x(node1);
        x2 = Body%x(node2);
        y1 = Body%y(node1);
        y2 = Body%y(node2);
#if NDIM==MDIM
        z1 =   Body%z(node1);
        z2 =   Body%z(node2);
#else
        z1=0.;
        z2=0.;
#endif
        vec(:)=(/x1-x2,y1-y2,z1-z2/);
        Body%Dist(i)= sqrt((vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)))*1.e-6  ! Real units
        
     end do
     ! Average distance 
     Dist_avg=sum(Body%Dist(:))/Body%nLinks;
     write(*,*) ' '
     write(*,*)'Dist_Avg',dist_avg, ' m'
     write(*,*) ' '
     Body%Lo_M(:)=Body%dist(:)/r_m;

  case (MATERIAL_FENE_POW , MATERIAL_WLC_POW , MATERIAL_WLC_FENE ,  MATERIAL_WLC_CQ) 

     body%lo=sqrt(((sm_rbc_ar*4.0)/(sqrt(3.)*body%nele)))*1.e-6;  ! in microns
     body%lo_m=body%lo/r_m;
  end select
  
  
  select case (sp_modeltype)
      
  case (MATERIAL_SF_FENE_POW)         
     
     !paramters in real units
     xo=1./body%r_mult;
     body%ks=(4 * mu_op)/(sqrt(3.)*( (2*xo*xo/(1-xo*xo)**2) +  &
          ((m_exp+1)/(1-xo*xo) )))
     
     !scaling the system;
     body%lmax_m=body%r_mult*body%lo_m;
     body%ks_m=body%ks*r_m/n_m
     
     do i=1,body%nlinks
        body%kp_m(i)=((body%ks(1)*body%dist(i)**real(m_exp+1))/(1-xo*xo))/(n_m * r_m * r_m)
    
     end do


#ifdef  PRINTOUT_PARMS
!!$     write(*,*)' '
!!$     write(*,*) 'lo= ',body%lo,' m'
!!$     write(*,*)' '
     write(*,*)'the fene spring coeff ks= ',body%ks,'n/m'
     write(*,*)' '
!!$     write(*,*) 'lo= ',body%lo_m,' dpd units'
     write(*,*)' '
     write(*,*) 'ks_m ', body%ks_m,' dpd units'
#endif
     
     
  case (MATERIAL_SF_WLC_POW)
     !******************************************************************
     !     
     !     worm-like chain + pow model 
     !
     !******************************************************************
     write(*,*)'----------------------------------------'
     write(*,*)'parameters for WLC_POW stress-free model'
     write(*,*)'----------------------------------------'
     
     !     paramters in real units

 
     xo=1./body%r_mult;
     
     body%p=(1./mu_op)*(sqrt(3.)*kboltz*sm_rbc_t/4.)*                        &
          (1./dist_avg)*(((xo/(2. *(1-xo)**3.))-(1./(4.*(1-xo)**2.))+0.25)+  &
          real(m_exp+1.)* ( ( 1. / ( 4. * (1-xo)**2. ) ) ) + xo - 0.25 ) 
     
     !     scaling the system;
     
     body%lmax_m=body%r_mult*body%lo_m;
     
     body%p_m=body%p/r_m;
        
     do i=1,body%nlinks
        body%kp_m(i)=((kboltz*sm_rbc_t*body%dist(i)**real(m_exp))/body%p)*  &
             ((0.25/((1.-xo)*(1.-xo))) + xo -0.25)/(n_m*r_m**real(m_exp))
     end do

#ifdef  PRINTOUT_PARMS
     write(*,*) 'rmult',body%r_mult
     write(*,*) ' '
     write(*,*) 'mu_op',mu_op
     write(*,*) ' '
     write(*,*) 'm_exp',m_exp
     write(*,*) ' '
     write(*,*) 'Dist_avg',dist_Avg
     write(*,*)' '
     write(*,*)'p= ',body%p,' m'
     write(*,*)' '  
     write(*,*) 'p_m ', body%p_m
#endif
     
     
     
!*******************************************************************

      case (MATERIAL_FENE_POW)

!*******************************************************************

!!$         body%lo=sqrt(((sm_rbc_ar*4.0)/(sqrt(3.)*body%nele)))*1e-6;   !in m 
!!$         write(12,*) 'lo= ',body%lo*1e6,' in micro-m'
!!$         
         !paramters in real units
        
         xo=1./2.05;
         
         body%ks_m=(4.*mu_op)/(sqrt(3.)*( (2.*xo*xo/(1.-xo*xo)**2.) + &
              ((m_exp+1.)/(1.-xo*xo) )))*r_m/n_m  ! n/m
         
         body%kp=(body%ks*body%lo**real(m_exp+1))/(1.-xo*xo)
         
         write(*,*)'ks   ',body%ks
         write(*,*)'kp_p ',body%kp

         
         !write(12,*)' ks= ',ks
         
         
         
         !scaling the system;
         !body%lo_m=body%lo/r_m;
         body%lmax_m=body%r_mult*body%lo_m
         
         body%ks_m=ks*r_m/n_m
         body%kp_m=(body%ks*body%dist**real(m_exp+1.))/(1.-xo*xo)/(n_m * r_m * r_m)           

        
      case (MATERIAL_WLC_CQ)

!******************************************************************
!     
!     worm-like chain + c model parameters
!
!******************************************************************
         
         
      case default
         call driver_abortflash (" no defined spring model was selected ")
         
      end select 
      


   !    element model types 
   sp_modeltype = body%mattype(1);
   
   select case (sp_modeltype)
      
   case (MATERIAL_SF_AREA_VOL ) 
      
      call gettri_norm(body%qms(:,1),body%tri,body%normals,body%centers,body%aoj,body%vot,nnodes,nele)
!!$      do i=1,nele
!!$         write(*,*) i,body%aoj(i)
!!$      end do
      !stop
      body%aot=sum(body%aoj(1:nele))
      body%at = body%aot;
      body%vt = body%vot;
      
   case(MATERIAL_AREA_VOL) 
      
   end select
   
   
   deallocate(normals, centers,vol_centers, new_comm) 
end subroutine sm_rbc_calc_parameters
	
