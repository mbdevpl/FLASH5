Module sm_pk_interface

  implicit none

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

  interface 
     subroutine sm_pk_init    
       implicit none
     end subroutine sm_pk_init
  end interface

  interface 
     subroutine sm_pk_finalize    
       implicit none
     end subroutine sm_pk_finalize
  end interface

  interface
      subroutine sm_pk_apply(ibd, time)
        implicit none
        integer, intent(in) :: ibd
        real,    intent(in) :: time
      end subroutine sm_pk_apply
  end interface

  interface
      subroutine sm_pk_ReadHDF5(infile_string)
        implicit none
        character, intent(in) :: infile_string*(*)
      end subroutine sm_pk_ReadHDF5
  end interface

  interface
      subroutine sm_pk_Flapping_BermanWang(time,n,Xe,v,vd,vdd,parameters)
        implicit none
        real, intent(in)    :: time
        integer, intent(in) :: n
        real, intent(in),  dimension(NDIM,n) :: Xe
        real, intent(in),  dimension(:)      :: parameters
        real, intent(out), dimension(NDIM,n) :: v, vd, vdd
      end subroutine sm_pk_Flapping_BermanWang
  end interface

  interface
      subroutine sm_pk_flapping_BermanWangFunctions(t, phi, theta, eta, parameters)
        implicit none
        real, intent(in)                :: t
        real, intent(out), dimension(NDIM) :: phi, theta, eta
        real, intent(in), dimension(:)  :: parameters
      end subroutine sm_pk_flapping_BermanWangFunctions
  end interface

  interface
      subroutine sm_pk_flapping_BermanWang_RotMat(phi, theta, eta, R, Rd, Rdd)
        implicit none
        real, intent(in),  dimension(MDIM)   :: phi, theta, eta
        real, intent(out), dimension(MDIM,MDIM) :: R, Rd, Rdd
      end subroutine sm_pk_flapping_BermanWang_RotMat
  end interface

  interface
      subroutine sm_pk_fixed(time,n,Xe,v,vd,vdd,parameters)
        implicit none
        real, intent(in)    :: time
        integer, intent(in) :: n
        real, intent(in),  dimension(NDIM,n) :: Xe
        real, intent(in),  dimension(:)      :: parameters
        real, intent(out), dimension(NDIM,n) :: v, vd, vdd
      end subroutine sm_pk_fixed
  end interface

  interface
     subroutine sm_pk_apply_entireBody(ibd, time)
       implicit none
       integer, intent(in) :: ibd
       real, intent(in) :: time
     end subroutine sm_pk_apply_entireBody
  end interface

  interface
     subroutine sm_pk_Whirl(time,n,Xe,v,vd,vdd,parameters)
        implicit none
        real, intent(in)    :: time
        integer, intent(in) :: n
        real, intent(in),  dimension(NDIM,n) :: Xe
        real, intent(in),  dimension(:)      :: parameters
        real, intent(out), dimension(NDIM,n) :: v, vd, vdd
      end subroutine sm_pk_Whirl
  end interface

  interface
     subroutine sm_pk_Whirl_functions(time, theta, parameters)
       implicit none
       real, intent(in)                     :: time
       real, intent(out), dimension(MDIM)   :: theta
       real, intent(in),  dimension(:)      :: parameters       
     end subroutine sm_pk_Whirl_functions
  end interface

  interface
     subroutine sm_pk_Whirl_RotMat(theta, nhat, R, Rd, Rdd)
       implicit none
       real, intent(in),  dimension(MDIM)      :: theta
       real, intent(in),  dimension(MDIM)      :: nhat
       real, intent(out), dimension(MDIM,MDIM) :: R, Rd, Rdd
     end subroutine sm_pk_Whirl_RotMat
  end interface

  interface
    subroutine sm_pk_fixed_dof(time,maxrestparams,paramcoord,vc,vcd,vcdd) 
      implicit none
      integer, intent(in) :: maxrestparams
      real, intent(in)    :: time, paramcoord(maxrestparams)
      real, intent(out)   :: vc, vcd, vcdd
    end subroutine sm_pk_fixed_dof
  end interface

  interface
     subroutine sm_pk_harmonic_dof(time,maxrestparams,paramcoord,vc,vcd,vcdd)
       implicit none
       integer, intent(in) :: maxrestparams
       real, intent(in)    :: time, paramcoord(maxrestparams)
       real, intent(out)   :: vc, vcd, vcdd
     end subroutine sm_pk_harmonic_dof
  end interface

  interface
     subroutine sm_pk_constvel_dof(time,maxrestparams,paramcoord,vc,vcd,vcdd)
       implicit none
       integer, intent(in) :: maxrestparams
       real, intent(in)    :: time, paramcoord(maxrestparams)
       real, intent(out)   :: vc, vcd, vcdd
     end subroutine sm_pk_constvel_dof
  end interface

  interface
     subroutine sm_pk_masterslave_rigid(maxdofs,nfix,xm,qmstr,qdmstr,qddmstr, &
                                        TNB,NwB_N,NaB_N,Xi,XiB,v,vd,vdd)
       implicit none
       integer, intent(in) :: nfix,maxdofs
       real, intent(in)  :: xm(MDIM,1)
       real, intent(in)  :: qmstr(maxdofs),qdmstr(maxdofs),qddmstr(maxdofs)
       real, intent(in)  :: TNB(MDIM,MDIM),NwB_N(MDIM,1),NaB_N(MDIM,1)
       real, intent(in)  :: Xi(MDIM,nfix),XiB(MDIM,nfix)
       real, intent(out) :: v(maxdofs,nfix),vd(maxdofs,nfix),vdd(maxdofs,nfix)
     end subroutine sm_pk_masterslave_rigid
  end interface

  interface
     subroutine sm_pk_Shaker(time,n,Xe,v,vd,vdd,parameters)
        implicit none
        real, intent(in)    :: time
        integer, intent(in) :: n
        real, intent(in),  dimension(NDIM,n) :: Xe
        real, intent(in),  dimension(:)      :: parameters
        real, intent(out), dimension(NDIM,n) :: v, vd, vdd
      end subroutine sm_pk_Shaker
  end interface

  interface
     subroutine sm_pk_Shaker_harmonic(time, omega, phi, Amp, tau, qe, qde, qdde)
       implicit none
       real, intent(in) :: time, omega, phi, Amp, tau
       real, intent(out) :: qe, qde, qdde
     end subroutine sm_pk_Shaker_harmonic
  end interface

  interface
     subroutine sm_pk_updatekinematics_rigid(ibd,restart,procedure_flag)
       implicit none 
       integer, intent(in) :: ibd,procedure_flag
       logical, intent(in) :: restart
     end subroutine sm_pk_updatekinematics_rigid
  end interface

  interface
     subroutine sm_pk_angvelconstraint_rigid(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_pk_angvelconstraint_rigid
  end interface


end Module sm_pk_interface
