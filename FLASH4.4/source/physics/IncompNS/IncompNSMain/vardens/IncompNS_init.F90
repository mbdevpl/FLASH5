!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS_init
!!
!! NAME
!!
!!  IncompNS_init
!!
!!
!! SYNOPSIS
!!
!!  call IncompNS_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine IncompNS_init(restart)

  use IncompNS_data
  use ImBound_interface, ONLY : ImBound_setData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm, Driver_getNstep

  use Multiphase_interface, ONLY : Multiphase_init 

  implicit none
  include 'Flash_mpi.h'
#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"
  logical, intent(IN) :: restart

  call Driver_getMype(MESH_COMM, ins_meshMe)
  call Driver_getNumProcs(MESH_COMM, ins_meshNumProcs)
  call Driver_getComm(MESH_COMM, ins_meshComm)

  call RuntimeParameters_get("cflflg", ins_cflflg)
  call RuntimeParameters_get("cfl", ins_cfl)
  call RuntimeParameters_get("isgs",ins_isgs)
  call RuntimeParameters_get("invRe",ins_invRe)
  call RuntimeParameters_get("sigma",ins_sigma)
  call RuntimeParameters_get("dtspec",ins_dtspec)
  call RuntimeParameters_get("intschm",ins_intschm)
  call RuntimeParameters_get("gravX",ins_gravX)
  call RuntimeParameters_get("gravY",ins_gravY)
  call RuntimeParameters_get("gravZ",ins_gravZ)
  call RuntimeParameters_get("dampC",ins_dampC)
  call RuntimeParameters_get("xDampL",ins_xDampL)
  call RuntimeParameters_get("zDampL",ins_zDampL)
  call RuntimeParameters_get("iConvU",ins_iConvU)
  call Driver_getNstep(ins_nstep)
  ins_restart=restart

  call RuntimeParameters_get("pressure_correct",ins_prescorr)

  !- kpd - This is altered to avoid pressure correction scheme
  if (ins_meshMe .eq. MASTER_PE) print*,"KPD - ins_prescorr is set to .FALSE. in IncompNS_init.F90"
  if (ins_meshMe .eq. MASTER_PE) print*,"         to avoid pressure correction scheme."
  ins_prescorr = .FALSE.

  ins_prescoeff = 0.
  if (ins_prescorr) ins_prescoeff = 1.

  call RuntimeParameters_get("vel_prolong_method",ins_prol_method)

  if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) 'ins_cfl   =',ins_cfl
     write(*,*) 'ins_isgs  =',ins_isgs
     write(*,*) 'ins_invRe =',ins_invRe
     write(*,*) 'ins_sigma =',ins_sigma
     write(*,*) 'ins_dtspec=',ins_dtspec
     write(*,*) 'ins_intschm=',ins_intschm
     write(*,*) 'ins_prescoeff=',ins_prescoeff
     write(*,*) 'vel_prolong_method=',ins_prol_method
     write(*,*) 'ins_gravX =',ins_gravX
     write(*,*) 'ins_gravY =',ins_gravY
     write(*,*) 'ins_gravZ =',ins_gravZ
     write(*,*) 'ins_dampC =',ins_dampC
     write(*,*) 'ins_iConvU =',ins_iConvU
  endif

  ! Call multiphase variables initialization routine:
  call Multiphase_init()	

  if(ins_meshMe==MASTER_PE)print*,'Incmp_Navier_Stokes initialized'

  ! Define value of kinematic viscosity for IB distributed forces:
  call ImBound_setData()

end subroutine IncompNS_init

