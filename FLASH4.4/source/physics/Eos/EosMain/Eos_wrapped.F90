!!****if* source/physics/Eos/EosMain/Eos_wrapped
!! NAME
!!
!!  Eos_wrapped
!! 
!! SYNOPSIS
!!
!!  call Eos_wrapped(  integer(IN) :: mode,
!!                     integer(IN) :: range(HIGH, MDIM),
!!                     integer(IN) :: blockID,
!!            optional,integer(IN) :: gridDataStruct )
!!
!! DESCRIPTION
!!
!! This function is provided for the user's convenience and acts as a simple
!! wrapper to the Eos interface. The Eos interface uses a single, flexible data
!! structure "eosData" to pass the thermodynamic quantities in and out of the
!! funtion (see Eos). The wrapper hides formation and use of eosData
!! from the users.
!!
!! While Eos does not know anything about blocks, Eos_wrapped takes its
!! input thermodynamic state variables from a given block's storage area.
!! It works by taking a selected section of a block described by array
!! "range" and translating it to eosData before calling the Eos routine.
!! Upon return from Eos, Eos_wrapper updates certain state variables in
!! the same section of the block's storage area. Which variables are taken
!! as input, and which are updated, depends on the "mode" argument.
!!
!! If you want to return the derived quantities defined from EOS_VAR+1:EOS_NUM
!! in Eos.h, then you must use the direct interface Eos().
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   mode : determines which variables are used as Eos input.
!!          Valid values are MODE_DENS_EI (where density and internal
!!          energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!          MODE_DENS_TEMP (density and temperature are inputs).
!!          These quantities are defined in constants.h, the argument is 
!!          forwarded unchanged to the Eos function call.
!!          Note that internal energy is grid variable EINT_VAR, not ENER_VAR.
!!
!! 
!!   range: an array that holds the lower and upper indices of the section
!!          of block on which Eos is to be applies. The example shows how
!!          the array describes the block section.
!!
!!   blockID: current block number
!!
!!   gridDataStruct : the grid data structure on whose data Eos is to be applied
!!
!!
!!  EXAMPLE 
!!      if range(LOW,IAXIS)=1,range(HIGH,IAXIS)=iguard,
!!         range(LOW,JAXIS)=1,range(HIGH,JAXIS)=jguard,
!!         range(LOW,KAXIS)=1,range(HIGH,KAXIS)=kguard,
!!      then Eos is applied to the lower left hand corner of the guard
!!      cells in the block. 
!!
!!      However, if the value were
!!         range(LOW,IAXIS)=iguard+1,range(HIGH,IAXIS)=iguard+nxb,
!!         range(LOW,JAXIS)=jguard+1,range(HIGH,JAXIS)=jguard+nyb,
!!         range(LOW,KAXIS)=kguard+1,range(HIGH,KAXIS)=kguard+nzb,
!!      then Eos is applied to all the interior cells in the block.
!!
!!  NOTES
!!      This interface is defined in Fortran Module 
!!      Eos_interface. All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_wrapped
!!
!!      This routine cannot use "INTERIOR" mode of indexing the range.  In the
!!      second example given above, although only the interior cells are being
!!      calculated with EOS, the range indices still must include the guard cells.
!!      See, for example, IsentropicVortex/Simulation_initBlock where the data is
!!      generated on INTERIOR cells with Grid_putRowData, but the same indices can't
!!      be used for the EOS call.
!!
!!  SEE ALSO
!!
!!     Eos
!!     Eos.h
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


!!$subroutine Eos_wrapped(mode,range,blockID, gridDataStruct)
!!$
!!$  use Driver_interface, ONLY : Driver_abortFlash
!!$  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
!!$  use Logfile_interface, ONLY: Logfile_stampMessage 
!!$  use Eos_interface, ONLY : Eos, Eos_putData, Eos_getData, Eos_arrayWrapped
!!$  use Eos_data, ONLY : eos_threadWithinBlock
!!$  !$ use omp_lib
!!$  implicit none
!!$
#include "Eos.h"
#include "constants.h"
#include "Flash.h"
!!$
!!$  integer, intent(in) :: mode
!!$  integer, dimension(2,MDIM), intent(in) :: range
!!$  integer,intent(in) :: blockID
!!$  integer, optional, intent(IN) :: gridDataStruct
!!$
!!$  real, pointer:: solnData(:,:,:,:)
!!$
!!$#ifndef FIXEDBLOCKSIZE
!!$  real, allocatable :: eosData(:),massFraction(:)
!!$#else
!!$  real, dimension(NSPECIES*MAXCELLS) :: massFraction
!!$  real, dimension(EOS_NUM*MAXCELLS) :: eosData
!!$#endif
!!$
!!$  logical,target,dimension(EOS_VARS+1:EOS_NUM) :: eosMask
!!$
!!$  integer :: ierr, dataStruct
!!$  integer :: i,j,k, vecLen
!!$  integer,dimension(MDIM) :: pos
!!$
!!$
!!$
!!$  call Grid_getBlkPtr(blockID,solnData,gridDataStruct)
!!$  call Eos_arrayWrapped(mode,range,solnData, gridDataStruct)
!!$
!!$  call Grid_releaseBlkPtr(blockID,solnData,gridDataStruct)
!!$
!!$  return
!!$end subroutine Eos_wrapped

!!****if* source/physics/Eos/EosMain/Eos_Wrapped
!! NAME
!!
!!  Eos_Wrapped
!! 
!! SYNOPSIS
!!
!!  call Eos_Wrapped(  integer(IN) :: mode,
!!                          integer(IN) :: range(HIGH, MDIM),
!!                     real,pointer(IN) :: solnData,
!!                 optional,integer(IN) :: gridDataStruct )
!!
!! DESCRIPTION
!!
!! This function is provided for the user's convenience and acts as a simple
!! wrapper to the Eos interface. The Eos interface uses a single, flexible data
!! structure "eosData" to pass the thermodynamic quantities in and out of the
!! funtion (see Eos). The wrapper hides formation and use of eosData
!! from the users.
!!
!! While Eos does not know anything about blocks, Eos_arrayWrapped takes its
!! input thermodynamic state variables from a given block's storage area.
!! It works by taking a selected section of a block
!! described by array "range" and translating it to eosData
!! before calling the Eos function.
!! Upon return from Eos, Eos_wrapper updates certain state variables in
!! the same section of the block's storage area. Which variables are taken
!! as input, and which are updated, depends on the "mode" argument.
!!
!! If you want to return the derived quantities defined from EOS_VAR+1:EOS_NUM
!! in Eos.h, then you must use the direct interface Eos().  Note that 
!! entropy EOS_ENTR is considered a derived variable.
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   mode : determines which variables are used as Eos input.
!!          The valid values are MODE_DENS_EI (where density and internal
!!          energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!          MODE_DENS_TEMP (density and temperature are inputs).
!!          These quantities are defined in constants.h, the argument is 
!!          forwarded unchanged to the Eos function call.
!!          Note that internal energy is grid variable EINT_VAR, not ENER_VAR.
!!
!! 
!!   range: an array that holds the lower and upper indices of the section
!!          of block on which Eos is to be applies. The example shows how
!!          the array describes the block section.
!!
!!   solnData: pointer to array section of solution data for a block
!!
!!   gridDataStruct : the grid data structure on whose data Eos is to be applied
!!
!!
!!  EXAMPLE 
!!      if range(LOW,IAXIS)=1,range(HIGH,IAXIS)=iguard,
!!         range(LOW,JAXIS)=1,range(HIGH,JAXIS)=jguard,
!!         range(LOW,KAXIS)=1,range(HIGH,KAXIS)=kguard,
!!      then Eos is applied to the lower left hand corner of the guard
!!      cells in the block. 
!!
!!      However if the value were
!!         range(LOW,IAXIS)=iguard+1,range(HIGH,IAXIS)=iguard+nxb,
!!         range(LOW,JAXIS)=jguard+1,range(HIGH,JAXIS)=jguard+nyb,
!!         range(LOW,KAXIS)=kguard+1,range(HIGH,KAXIS)=kguard+nzb,
!!      then Eos is applied to all the interior cells in the block.
!!
!!  NOTES
!!      This interface is defined in Fortran Module 
!!      Eos_interface. All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_arrayWrapped
!!
!!      This routine cannot use "INTERIOR" mode of indexing the range.  In the
!!      second example given above, although only the interior cells are being
!!      calculated with EOS, the range indices still must include the guard cells.
!!      See, for example, IsentropicVortex/Simulation_initBlock where the data is
!!      generated on INTERIOR cells with Grid_putRowData, but the same indices can't
!!      be used for the EOS call.
!!
!!  SEE ALSO
!!
!!     Eos
!!     Eos.h
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Eos_wrapped(mode,range,solnData, gridDataStruct)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Logfile_interface, ONLY: Logfile_stampMessage 
  use Eos_interface, ONLY : Eos, Eos_putData, Eos_getData
  use Eos_data, ONLY : eos_threadWithinBlock
  implicit none

#include "FortranLangFeatures.fh"

  integer, intent(in) :: mode
  integer, dimension(2,MDIM), intent(in) :: range
  real, POINTER_INTENT_IN :: solnData(:,:,:,:)
  integer, optional, intent(IN) :: gridDataStruct


  real, allocatable :: eosData(:),massFraction(:)

  logical,target,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: ierr, dataStruct
  integer :: i,j,k, vLen1D, vecLen,offset
  integer,dimension(MDIM) :: pos


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
     call Driver_abortFlash("[Eos_arrayWrapped] "//&
          "invalid mode: must be MODE_DENS_PRES, MODE_DENS_TEMP, MODE_DENS_EI, or variants thereof, or MODE_EOS_NOP")
  end if
#endif

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *

  vLen1D = range(HIGH,IAXIS)-range(LOW,IAXIS)+1
  vecLen=vLen1D*(range(HIGH,JAXIS)-range(LOW,JAXIS)+1)*(range(HIGH,KAXIS)-range(LOW,KAXIS)+1)
  if (vecLen==0) return ! * Return immediately for empty IAXIS range! (for efficiency and avoiding index range errors)

  ! solnData points to solution data in UNK (or other data structure).
  ! The length of the data being operated upon is determined from the range input array.

  if(present(gridDataStruct))then
     dataStruct=gridDataStruct
  else
     dataStruct=CENTER
  end if

  allocate(massFraction(NSPECIES*vecLen))
  allocate(eosData(EOS_NUM*vecLen))

  eosMask = .FALSE.

  pos(IAXIS)=range(LOW,IAXIS)

  call Eos_getData(range,vecLen,solnData,dataStruct,eosData,massFraction)
  call Eos(mode,vecLen,eosData,massFraction,mask=eosMask)

  call Eos_putData(range,vecLen,solnData,dataStruct,eosData)
        
  deallocate(eosData)
  deallocate(massFraction)

  return
end subroutine Eos_wrapped
