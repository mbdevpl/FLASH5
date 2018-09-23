!!****if* source/Simulation/SimulationMain/unitTest/Eos/Driver_evolveFlash
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
!!  This is the main for the Eos unit test
!!
!!***

subroutine Driver_evolveFlash()

#include "constants.h"

  use Driver_data, ONLY: dr_globalMe 
  use Eos_interface, ONLY : Eos_unitTest
  use Logfile_interface,   ONLY : Logfile_stamp, Logfile_close
  use Timers_interface,    ONLY : Timers_getSummary
  use IO_interface,        ONLY : IO_writeCheckpoint, IO_outputFinal
  use Grid_interface, ONLY : Grid_getMaxRefinement, Grid_getBlkPtr,Grid_releaseBlkPtr
  use gr_interface, ONLY : gr_getBlkIterator, gr_releaseBlkIterator
  use gr_iterator, ONLY : gr_iterator_t
  use block_metadata, ONLY : block_metadata_t

  implicit none

  interface
     integer function ut_getFreeFileUnit()
     end function ut_getFreeFileUnit
  end interface

  logical :: perfect,thisBlock

  character(len=20) :: fileName
  integer           :: fileUnit
  integer,dimension(4) :: prNum
  integer :: temp,i
  integer,dimension(LOW:HIGH,MDIM)::tileLimits
  type(gr_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc
  integer:: level, maxLev
  real,pointer,dimension(:,:,:,:) :: Uout


  ! stays true if no errors are found
  perfect = .true.
  
  call Logfile_stamp( 'Entering testing loop' , '[Driver_evolveFlash]')

  temp = dr_globalMe

  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))

  fileUnit = ut_getFreeFileUnit()
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe

  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
  do level=1,maxLev
     call gr_getBlkIterator(itor, level=level)
     do while(itor%is_valid())
        call itor%blkMetaData(blockDesc)
        
        tileLimits = blockDesc%limits
        call Grid_getBlkPtr(blockDesc, Uout)

        print *, "Preparing to call Eos_unitTest"
        thisBlock=.true.
        Call Eos_unitTest(fileUnit,thisBlock, Uout, tileLimits, blockDesc)
        print *, "    Returned from Eos_unitTest"
        perfect=perfect.and.thisBlock
        call Grid_releaseBlkPtr(blockDesc, Uout)
        call itor%next()
     end do
     call gr_releaseBlkIterator(itor)
  end do
  
  if (perfect) then
     write(fileUnit,'("all results conformed with expected values.")')
  endif
  
  close(fileUnit)
  
  call Logfile_stamp( 'Exiting testing loop' , '[Driver_evolveFlash]')
  !if we write files, do final output
  call IO_writeCheckpoint()
  call IO_outputFinal()
  call Timers_getSummary( 0)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()
  
  return
  
end subroutine Driver_evolveFlash
      
