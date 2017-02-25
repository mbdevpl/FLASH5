!!****if* source/physics/Eos/EosMain/Eos_everywhere
!! NAME
!!
!!  Eos_everywhere
!!
!! SYNOPSIS
!!
!!  call Eos_everywhere(  integer(IN) :: mode,
!!               optional,integer(IN) :: gridDataStruct )
!!
!! DESCRIPTION
!!
!! Apply Eos() to all blocks.
!!
!!  ARGUMENTS
!!
!!   mode : determines which variables are used as Eos input.
!!          The valid values are MODE_DENS_EI (where density and internal
!!          energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!          MODE_DENS_TEMP (density and temperature are inputs).
!!          These quantities are defined in constants.h, the argument is
!!          forwarded unchanged to the Eos function call.
!!
!!   gridDataStruct : the grid data structure on whose data Eos is to be applied
!!
!!
!!  EXAMPLE
!!
!!  NOTES
!!
!!  SEE ALSO
!!
!!     Eos
!!     Eos_wrapped
!!     Grid_makeVector
!!
!!***

#include "constants.h"
#include "Eos.h"


subroutine Eos_everywhere(mode,gridDataStruct)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : GRID_COPYDIR_TO_VECT, GRID_COPYDIR_FROM_VECT
  use Grid_interface,   ONLY : Grid_makeVector
  use Eos_interface,    ONLY : Eos
  implicit none

  integer, intent(in) :: mode
  integer, optional, intent(IN) :: gridDataStruct

  integer,parameter :: VECTOR_LEN_TEST = 80
  integer,parameter :: NUM_VECTORS_TEST = 60000

  real, allocatable :: eosVectors(:,:,:)
!!$  real,target :: eosVectors(VECTOR_LEN_TEST,EOS_NUM,NUM_VECTORS_TEST)

  logical,target,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: ierr, dataStruct
  integer :: i, vecLen,vecNum, vecLastFree


!! ---------------------------------------------------------------------------------
  ! Test calling arguments
#define DEBUG
#ifdef DEBUG
  ierr = 1
  select case (mode)
  case (MODE_DENS_PRES)
     ierr = 0
  case (MODE_DENS_TEMP)
     ierr = 0
  case (MODE_DENS_EI)
     ierr = 0
  case (MODE_EOS_NOP)
     ierr = 0
  case (MODE_DENS_TEMP_ALL,MODE_DENS_TEMP_EQUI)
     ierr = 0
  case (MODE_DENS_EI_ALL,MODE_DENS_EI_SCATTER,MODE_DENS_EI_GATHER)
     ierr = 0
  case (MODE_DENS_EI_SELE_GATHER)
     ierr = 0
  case (MODE_DENS_ENTR)
     ierr = 0
  end select

  if(ierr /= 0) then
     call Driver_abortFlash("[Eos_everywhere] "//&
          "invalid mode: must be MODE_DENS_PRES, MODE_DENS_TEMP, MODE_DENS_EI, or variants thereof, or MODE_EOS_NOP")
  end if
#endif

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *


  if(present(gridDataStruct))then
     dataStruct=gridDataStruct
  else
     dataStruct=CENTER
  end if

  vecLen = VECTOR_LEN_TEST
  vecNum = NUM_VECTORS_TEST
  allocate(eosVectors(vecLen,EOS_NUM,vecNum))

  eosMask = .FALSE.

  call Grid_makeVector(vecLen,EOS_NUM,eosVectors,vecNum,vecLastFree, copyDirection=GRID_COPYDIR_TO_VECT,&
       gridDataStruct=dataStruct)
  if (vecLastFree < 0) then
     call Driver_abortFlash("Eos_everywhere: number of vectors for Grid_makeVector too small!")
  end if


  do i=1, vecNum
     if (i < vecNum .OR. vecLastFree == 0) then
        call Eos(mode,vecLen,eosVectors(1,1,i),mask=eosMask)
     else
        call Eos(mode,vecLen,eosVectors(1,1,i),mask=eosMask, vecEnd=VECTOR_LEN_TEST-veclastFree)
     end if
  end do

  call Grid_makeVector(vecLen,EOS_NUM,eosVectors,vecNum, copyDirection=GRID_COPYDIR_FROM_VECT,&
       gridDataStruct=dataStruct)

  deallocate(eosVectors)


  return
end subroutine Eos_everywhere
