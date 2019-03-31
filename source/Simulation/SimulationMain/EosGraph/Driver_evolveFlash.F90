!!****if* source/Simulation/SimulationMain/EosGraph/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  call Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! 
!! Call Eos for various temperatures and write results to files.
!!
!!  
!!***



#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

  use Eos_interface
#include "Eos.h"
#include "constants.h"
#include "Flash.h"
  implicit none
  integer :: i,j,k, outLUnitE, outLUnitP
  real,parameter :: xOgMin = 1e4
  real,parameter :: xOgMax = 1e11
  integer,parameter :: nIOg = 70
#ifdef FLASH_3T
  integer,parameter :: eosMode = MODE_DENS_TEMP_EQUI
#else
  integer,parameter :: eosMode = MODE_DENS_TEMP
#endif
  real :: eosData(EOS_NUM)
  logical :: eosMask(EOS_VARS+1:EOS_NUM)
  real,allocatable :: xOgFaces(:)
  real :: xOgStep, xOgFact
  real :: density,temp,x,y,y1,y2,y3

  eosMask(:) = .FALSE.
  eosMask(EOS_DET) = .TRUE.
  eosMask(EOS_DPT) = .TRUE.
#ifdef FLASH_3T
  eosMask(EOS_EINTION) = .TRUE.
  eosMask(EOS_EINTELE) = .TRUE.
  eosMask(EOS_EINTRAD) = .TRUE.
  eosMask(EOS_PRESION) = .TRUE.
  eosMask(EOS_PRESELE) = .TRUE.
  eosMask(EOS_PRESRAD) = .TRUE.
#else
  eosMask(EOS_DED) = .TRUE.
  eosMask(EOS_DPD) = .TRUE.
#endif

  allocate(xOgFaces(0:nIOg+2))

  if (xOgMin .LE. 0.0) then     !linear spacing...

     xOgStep = (xOgMax - xOgMin) / nIOg
     do i=0,nIOg+2
        xOgFaces(i) = xOgMin + (real(i)-1.0) * xOgStep 
     end do

  else                          !log spacing...
     xOgFact = 10.0**((alog10(xOgMax) - alog10(xOgMin)) / nIOg)
     do i=0,nIOg+2
        xOgFaces(i) = xOgMin * xOgFact**(real(i)-1.0)
     end do
  end if
  
  outLUnitE = 30
  open(outLUnitE,file='EOSdumpE.dat',status='UNKNOWN')
  outLUnitP = 31
  open(outLUnitP,file='EOSdumpP.dat',status='UNKNOWN')

  density = 1.0e3

  do i = 1,nIOg+1

     x = xOgFaces(i)
     eosData(EOS_DENS) = density
     eosData(EOS_TEMP) = x
     eosData(EOS_ABAR) = 1.0
     eosData(EOS_ZBAR) = 1.0
     call Eos(eosMode,1,eosData,mask=eosMask)

     y = eosData(EOS_EINT)
#ifdef FLASH_3T
     y1 = eosData(EOS_EINTION)
     y2 = eosData(EOS_EINTELE)
     y3 = eosData(EOS_EINTRAD)
#else
     y1 = eosData(EOS_DED)
     y2 = eosData(EOS_ENTR)
     y3 = eosData(EOS_ABAR)
#endif
     write(outLUnitE,1) x,y,y1,y2,y3,eosData(EOS_DET)

     y = eosData(EOS_PRES)
#ifdef FLASH_3T
     y1 = eosData(EOS_PRESION)
     y2 = eosData(EOS_PRESELE)
     y3 = eosData(EOS_PRESRAD)
#else
     y1 = eosData(EOS_DPD)
     y2 = eosData(EOS_GAMC)
     y3 = eosData(EOS_ZBAR)
#endif
     write(outLUnitP,1) x,y,y1,y2,y3,eosData(EOS_DPT)
1    format(6(1x,1PG12.3))

  end do

  return
  
end subroutine Driver_evolveFlash



