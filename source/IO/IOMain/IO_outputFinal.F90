!!****if* source/IO/IOMain/IO_outputFinal
!!
!! NAME
!!
!!  IO_outputFinal
!!
!!
!! SYNOPSIS
!!
!!  IO_outputFinal()
!!
!!
!! DESCRIPTION
!!
!!  This routine is called after the code has exited the main timestep
!!  loop.  It outputs the last checkpoint,plotfile and particle plotfile.
!!
!!  If particles are not included a stub (empty) routine will be called.
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  The state of module level logical variable io_outputInStack.
!!
!!***

subroutine IO_outputFinal()

#include "constants.h"  
  use IO_data, ONLY : io_justCheckpointed, io_outputInStack, &
       io_memoryStatFreq,  &
       io_summaryOutputOnly
  use IO_interface, ONLY : IO_writeCheckpoint, IO_writePlotfile, &
    IO_writeParticles

  implicit none


  !This setting is used to ensure valid data throughout grid and ancestor blocks
  io_outputInStack = .true.

  if (.not.io_summaryOutputOnly) then
     if(.not. io_justCheckpointed) then  
        call IO_writeCheckpoint()
     end if
   
     call IO_writePlotfile(.true.)
     
!! Devnote :: preprocessors because amrex particles are being handled through amrex grid
#ifndef FLASH_GRID_AMREX
     call IO_writeParticles(.false.)
#endif
  end if

  io_outputInStack = .false.

  !------------------------------------------------------------------------------
  ! Dump out memory usage statistics if we are monitoring them,
  ! after writing binary output files for the last time.
  !------------------------------------------------------------------------------
  if (io_memoryStatFreq > 0) call io_memoryReport()

end subroutine IO_outputFinal
