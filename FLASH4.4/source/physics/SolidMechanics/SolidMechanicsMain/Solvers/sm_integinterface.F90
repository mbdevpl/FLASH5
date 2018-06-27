

module sm_integinterface

#include "Flash.h"

  ! Wrapper Functions
  interface
     subroutine sm_integ_advance(restart_local)
       implicit none
       logical, INTENT(IN) :: restart_local
     end subroutine sm_integ_advance
  end interface

  interface
     subroutine sm_integ_checkconverg(convflag_all)
       implicit none
       integer, INTENT(OUT) :: convflag_all
     end subroutine sm_integ_checkconverg
  end interface

  interface
     subroutine sm_integ_adv1dt(restart_local)
       implicit none
       logical, INTENT(IN) :: restart_local
     end subroutine sm_integ_adv1dt
  end interface

  interface
     subroutine sm_integ_writeCheckpoint()
       implicit none
     end subroutine sm_integ_writeCheckpoint
  end interface

  ! Predictor-Corrector methods routines
  interface
      subroutine sm_PredCorr_init(ibd,restart)
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
      end subroutine sm_PredCorr_init
  end interface

  interface
      subroutine sm_PredCorr_advance(ibd,restart)
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
      end subroutine sm_PredCorr_advance
  end interface

  interface

     subroutine sm_predcorrInteg(ibd,dt)
       implicit none
       integer, intent(in) :: ibd
       real, intent(in) :: dt
     end subroutine sm_predcorrInteg

  end interface

  interface
     subroutine sm_PredCorr_checkconverg(testflg,ibd,convflag_pc,errmax_pc)
       implicit none
       integer, intent(IN)  :: testflg,ibd
       integer, intent(OUT) :: convflag_pc
       real, intent(OUT)    :: errmax_pc
     end subroutine sm_PredCorr_checkconverg
  end interface


  interface
     subroutine sm_pcem(n,X0,DX0,dt,X1)
       implicit none
       integer, intent(in) :: n
       real, intent(in) :: dt
       real, dimension(n), intent(in) :: X0,DX0
       real, dimension(n), intent(out):: X1
     end subroutine sm_pcem
  end interface

  interface

     subroutine sm_pcmem(n,X0,DX0,DX1,dt,X1)
       implicit none
       integer, intent(in) :: n
       real, intent(in) :: dt
       real, dimension(n), intent(in) :: X0,DX0,DX1
       real, dimension(n), intent(out):: X1
     end subroutine sm_pcmem

  end interface

  interface
     subroutine sm_pc_AB2(n,X0,DX0,DXm1,dt,m,X1)
       implicit none
       integer, intent(in) :: n, m
       real,    intent(in) :: dt(-m:0)
       real, dimension(n), intent(in) :: X0,DX0,DXm1
       real, dimension(n), intent(out):: X1
     end subroutine sm_pc_AB2
  end interface

  interface
     subroutine sm_pc_AM2(n,X0,DX0,DXm1,DX1,dt,m,X1)
       implicit none
       integer, intent(in) :: n, m
       real,    intent(in) :: dt(-m:0)
       real, dimension(n), intent(in) :: DXm1,X0,DX0,DX1
       real, dimension(n), intent(out):: X1
     end subroutine sm_pc_AM2
  end interface

  interface
     subroutine sm_pc_AB3(n,X0,DX0,DXm1,DXm2,dt,m,X1)
       implicit none
       integer, intent(in) :: n, m
       real,    intent(in) :: dt(-m:0)
       real, dimension(n), intent(in) :: X0,DX0,DXm1,DXm2
       real, dimension(n), intent(out):: X1
     end subroutine sm_pc_AB3
  end interface

  interface
     subroutine sm_pc_AM3(n,X0,DX0,DXm1,DXm2,DX1,dt,m,X1)
       implicit none
       integer, intent(in) :: n, m
       real,    intent(in) :: dt(-m:0)
       real, dimension(n), intent(in) :: DXm1,DXm2,X0,DX0,DX1
       real, dimension(n), intent(out):: X1
     end subroutine sm_pc_AM3
  end interface

  interface
     subroutine sm_pc_AB4(n,X0,DX0,DXm1,DXm2,DXm3,dt,m,X1)
       implicit none
       integer, intent(in) :: n, m
       real,    intent(in) :: dt(-m:0)
       real, dimension(n), intent(in) :: X0,DX0,DXm1,DXm2,DXm3
       real, dimension(n), intent(out):: X1
     end subroutine sm_pc_AB4
  end interface

  interface
     subroutine sm_pc_AM4(n,X0,DX0,DXm1,DXm2,DXm3,DX1,dt,m,X1)
       implicit none
       integer, intent(in) :: n, m
       real,    intent(in) :: dt(-m:0)
       real, dimension(n), intent(in) :: DXm1,DXm2,DXm3,X0,DX0,DX1
       real, dimension(n), intent(out):: X1
     end subroutine sm_pc_AM4
  end interface

  interface 
      subroutine sm_pc_initmass(ibd)
        implicit none
        integer, intent(in) :: ibd
      end subroutine sm_pc_initmass
  end interface

  interface
      subroutine sm_pc_initdamp(ibd)
        implicit none
        integer, intent(in) :: ibd
      end subroutine sm_pc_initdamp
  end interface

  interface
      subroutine sm_pc_compute_qddn(ibd)
        implicit none
        integer, intent(in) :: ibd
      end subroutine sm_pc_compute_qddn
  end interface

   interface
     subroutine sm_PredCorr_readCheckpoint(ibd)
       integer, intent(in) :: ibd
     end subroutine sm_PredCorr_readCheckpoint
  end interface

  interface
     subroutine sm_PredCorr_writeCheckpoint(ibd, sm_checkpt_num)
       integer, intent(in) :: ibd
       integer, intent(in) :: sm_checkpt_num
     end subroutine sm_PredCorr_writeCheckpoint
  end interface

  interface
     subroutine sm_pc_isPredictor(ibd, flag)
       implicit none
       integer, intent(in)  :: ibd
       integer, intent(out) :: flag
     end subroutine sm_pc_isPredictor
  end interface

  interface
     subroutine sm_pc_isCorrector(ibd, flag)
       implicit none
       integer, intent(in)  :: ibd
       integer, intent(out) :: flag
     end subroutine sm_pc_isCorrector
  end interface
  
  ! Verlet Integrator routines     
  interface
      subroutine sm_Verlet_init(ibd,restart)
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
      end subroutine sm_Verlet_init
  end interface

  interface
      subroutine sm_Verlet_advance(ibd,restart)
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
      end subroutine sm_Verlet_advance
  end interface

    interface 
     subroutine sm_verlet_partone(ibd,dt)       
       implicit none
       integer :: ibd
       real, intent(IN) :: dt
     end subroutine sm_verlet_partone
  end interface
  
  interface
     subroutine sm_verlet_parttwo(ibd,dt)    
       implicit none       
       integer :: ibd
       real, intent(IN) :: dt
     end subroutine sm_verlet_parttwo
  end interface

  interface
     subroutine sm_Verlet_writeCheckpoint(ibd)
       integer, intent(in) :: ibd
     end subroutine sm_Verlet_writeCheckpoint
  end interface

  !
  ! DT Checks:
  !
  interface
     subroutine sm_EstDT_3DFlexible(ibd, dt)
       implicit none
       integer, intent(in) :: ibd
       real, intent(out) :: dt
     end subroutine sm_EstDT_3DFlexible
  end interface

  interface
     subroutine sm_EstDT_rigid(ibd, dt)
       implicit none
       integer, intent(in) :: ibd
       real, intent(out) :: dt
     end subroutine sm_EstDT_rigid
  end interface

  interface
     subroutine sm_EstDT_rbc(ibd, dt)
       implicit none
       integer, intent(in) :: ibd
       real, intent(out) :: dt
     end subroutine sm_EstDT_rbc
  end interface

  ! Generalized Alpha Predictor-Corrector methods routines
  interface
      subroutine sm_GenAlpha_init(ibd,restart)
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
    end subroutine sm_GenAlpha_init
  end interface

  interface
      subroutine sm_GenAlpha_advance(ibd,restart)
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
    end subroutine sm_GenAlpha_advance
  end interface

  interface
     subroutine sm_GenAlpha_checkconverg(testflg,ibd,convflag_pc,errmax_pc)
       implicit none
       integer, intent(IN)  :: testflg,ibd
       integer, intent(OUT) :: convflag_pc
       real, intent(OUT)    :: errmax_pc
     end subroutine sm_GenAlpha_checkconverg
  end interface

  interface
     subroutine sm_ga_compute_qddn0(ibd)
       implicit none
       integer, intent(in) :: ibd ! body number
     end subroutine sm_ga_compute_qddn0
  end interface

  interface
     subroutine sm_ga_initdamp(ibd)
       implicit none
       integer, intent(in) :: ibd ! body number
     end subroutine sm_ga_initdamp
  end interface  

  interface
     subroutine sm_GenAlpha_readCheckpoint(ibd)
       integer, intent(in) :: ibd
     end subroutine sm_GenAlpha_readCheckpoint
  end interface

  interface
     subroutine sm_GenAlpha_writeCheckpoint(ibd, sm_checkpt_num )
       integer, intent(in) :: ibd
       integer, intent(in) :: sm_checkpt_num
     end subroutine sm_GenAlpha_writeCheckpoint
  end interface

  interface
     subroutine sm_ga_modstiff_fluid(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_ga_modstiff_fluid
  end interface

  interface
     subroutine sm_ga_isPredictor(ibd, flag)
       implicit none
       integer, intent(in)  :: ibd
       integer, intent(out) :: flag
     end subroutine sm_ga_isPredictor
  end interface

    interface
     subroutine sm_ga_isCorrector(ibd, flag)
       implicit none
       integer, intent(in)  :: ibd
       integer, intent(out) :: flag
     end subroutine sm_ga_isCorrector
  end interface


end module sm_integinterface
