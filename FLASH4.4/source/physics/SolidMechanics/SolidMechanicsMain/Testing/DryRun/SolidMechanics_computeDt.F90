!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Testing/DryRun/SolidMechanics_computeDt
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!  
!!
!!***

subroutine SolidMechanics_computeDt(dt_solid)
  use SolidMechanics_data, only: sm_BodyInfo
  use sm_integinterface, only: sm_EstDT_3DFlexible
  use Driver_interface, only: Driver_getNStep, Driver_getDT
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  real, INTENT(INOUT) :: dt_solid
  real :: dt

  integer :: nstep, ibd

  ibd = 1

  !dt_solid = 0.000001 !!! UBER-HACK !!  

  call Driver_getNStep(nstep)
  call Driver_getDT(dt)

  if( nstep > sm_BodyInfo(ibd)%ndofs ) then

     if( dt+1.e-5 <= 0.5e-4 ) then
        dt_solid = dt+1.e-5
     else
           dt_solid = dt
     end if

  else
     dt_solid = 1.E-5
  end if

!!$
!!$  if( nstep < 2*sm_BodyInfo(ibd)%ndofs ) then
!!$     dt_solid = 1.e-6
!!$     return
!!$  end if


  !if( mod(nstep,200) == 0 ) then
  !   call sm_EstDT_3DFlexible(ibd,dt)
  !   write(*,*) 'Est DT = ',dt
  !end if

  call RuntimeParameters_get("dtspec",dt_solid)
  !dt_solid = 0.002
  
  return

end subroutine SolidMechanics_computeDt
