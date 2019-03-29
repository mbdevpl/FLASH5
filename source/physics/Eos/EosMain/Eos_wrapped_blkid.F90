!!****if* source/physics/Eos/EosMain/Eos_wrapped_blkid
!! NAME
!!
!!  Eos_wrapped_blkid
!! 
!! SYNOPSIS
!!
!!  call Eos_wrapped_blkid(  integer(IN) :: mode,
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
!! While Eos does not know anything about blocks, Eos_wrapped_blkid takes its
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
!!      use Eos_interface, ONLY : Eos_wrapped_blkid
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


subroutine Eos_wrapped_blkid(mode,range,blockID, gridDataStruct)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Logfile_interface, ONLY: Logfile_stampMessage 
  use Eos_interface, ONLY : Eos, Eos_putData, Eos_getData, Eos_wrapped
  use Eos_data, ONLY : eos_threadWithinBlock
  !$ use omp_lib
  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: mode
  integer, dimension(2,MDIM), intent(in) :: range
  integer,intent(in) :: blockID
  integer, optional, intent(IN) :: gridDataStruct

  real, pointer:: solnData(:,:,:,:)

#ifndef FIXEDBLOCKSIZE
  real, allocatable :: eosData(:),massFraction(:)
#else
  real, dimension(NSPECIES*MAXCELLS) :: massFraction
  real, dimension(EOS_NUM*MAXCELLS) :: eosData
#endif

  logical,target,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: ierr, dataStruct
  integer :: i,j,k, vecLen
  integer,dimension(MDIM) :: pos



  call Grid_getBlkPtr(blockID,solnData,gridDataStruct)
  call Eos_wrapped(mode,range,solnData, gridDataStruct)

  call Grid_releaseBlkPtr(blockID,solnData,gridDataStruct)

  return
end subroutine Eos_wrapped_blkid
  
