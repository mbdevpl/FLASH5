!!****if* source/Simulation/SimulationMain/Sedov/Simulation_finalize
!!
!! NAME
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This dummy function cleans up the Simulation unit, deallocates memory, etc.
!!  However, as nothing needs to be done, only this stub is included.
!!
!! ARGUMENTS
!!
!!
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Simulation_finalize()
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_dump
  use Simulation_data, ONLY: sim_fileUnitOutNum, sim_fileUnitOutAna
  use Simulation_data, ONLY : sim_globalMe, &
        nLargestMaxSummary,  &
        sim_testInitialized, &
        sim_testFirstVals, &
        sim_testLastVals, &
        sim_testLargestVals, &
        sim_testLargestWhen

  implicit none


  character (len = MAX_STRING_LENGTH) :: fileName

  logical :: perfect, imperfect

  integer :: i
  integer :: fileUnit
  integer :: ut_getFreeFileUnit
  real    :: tols(4)

#if NDIM==1
  integer :: blkcnt
  integer :: blklst(MAXBLOCKS)

  integer :: lb

  call Grid_getListOfBlocks(LEAF,blklst,blkcnt)
  do lb = 1, blkcnt
     call Grid_dump((/DENS_VAR,PRES_VAR,VELX_VAR,EINT_VAR/),4, blklst(lb),gcell=.FALSE.)
  end do
#endif

  close(sim_fileUnitOutNum)
  close(sim_fileUnitOutAna)

  if (sim_globalMe .NE. MASTER_PE) return ! avoid writing lots of (empty) files if using many procs...

  perfect = .TRUE.

  tols(1) = 1.e-10
  tols(2) = 1.e-6
  tols(3) = 1.e-6
  tols(4) = 1.e-6

  fileUnit = ut_getFreeFileUnit ()
  fileName = "LargestSummary.out"

  open (fileUnit, file = fileName)
!!$  if (sim_globalMe==MASTER_PE) then
!!$     do i=1,min(4,1+N_DIM)
!!$        imperfect = (abs(sim_testLargestVals(i)-sim_testRefVals(i)) > tols(i))
!!$        perfect = (perfect .AND. .NOT. imperfect)
!!$        if (.NOT. perfect) then
!!$           select case (i)
!!$           case(1)
!!$              call Logfile_stamp(sim_testLargestWhen(i), &
!!$                   '[Simulation_finalize] Total mass conservation - worst at')
!!$           case(2)
!!$              call Logfile_stamp(sim_testLargestWhen(i), &
!!$                   '[Simulation_finalize] Total X momentum conservation - worst at')
!!$           case(3)
!!$              if (imperfect .OR. N_DIM>1) call Logfile_stamp(sim_testLargestWhen(i), &
!!$                   '[Simulation_finalize] Total Y momentum conservation - worst at')
!!$           case(4)
!!$              if (imperfect .OR. N_DIM>2) call Logfile_stamp(sim_testLargestWhen(i), &
!!$                   '[Simulation_finalize] Total Z momentum conservation - worst at')
!!$           end select
!!$        end if
!!$        if (imperfect) then
!!$           call Logfile_stamp(sim_testLargestVals(i)-sim_testRefVals(i), &
!!$             '[Simulation_finalize] Test failed,    Error')
!!$           call Logfile_stamp(sim_testFirstVals(i), &
!!$             '[Simulation_finalize] Test failed, reference value')
!!$           call Logfile_stamp(sim_testLastVals(i), &
!!$             '[Simulation_finalize] Test failed, last value')
!!$           call Logfile_stamp(sim_testLargestVals(i), &
!!$             '[Simulation_finalize] Test failed, worst value')
!!$        else if (.NOT. perfect .AND. i .LE. 1+N_DIM) then
!!$           call Logfile_stamp(sim_testLargestVals(i)-sim_testRefVals(i), &
!!$             '[Simulation_finalize] Test succeeded, Error')
!!$           call Logfile_stamp(sim_testFirstVals(i), &
!!$             '[Simulation_finalize] Test succeeded, reference value')
!!$           call Logfile_stamp(sim_testLastVals(i), &
!!$             '[Simulation_finalize] Test succeeded, last value')
!!$           call Logfile_stamp(sim_testLargestVals(i), &
!!$             '[Simulation_finalize] Test succeeded, worst value')
!!$        end if
!!$     end do
!!$  end if


!!$  if (perfect) then
!!$      write (fileUnit,'(a)') 'SUCCESS all results conformed with expected values.'
!!$  else
!!$      write (fileUnit,'(a)') 'FAILURE'
!!$  end if

!!$97 format(1x,"=========",1x,I3,1x,"=========")
97 format(1x,1x,I3,".  ")
98 format(1x,A,":",(2x,1PG14.7,6x))
99 format(1x,A,(1PG17.10))
  if (sim_globalMe==MASTER_PE) then
     do i=0,nLargestMaxSummary
        write (fileUnit,97,ADVANCE="NO") i+1
        write (fileUnit,98,ADVANCE="NO") "First val  ",sim_testFirstVals(i)
        write (fileUnit,98,ADVANCE="NO") "Last val   ",sim_testLastVals(i)
        write (fileUnit,98,ADVANCE="NO") "Largest val",sim_testLargestVals(i)
        write (fileUnit,99) "Largest at t=",sim_testLargestWhen(i)
     end do
  end if
  close (fileUnit)


  return

end subroutine Simulation_finalize
