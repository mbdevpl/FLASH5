!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/rigid/sm_ioInit_rigid
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
!!  Initialize IO variables and files for SolidMechanics rigid bodies.
!!
!!***

subroutine sm_ioInit_rigid(restart,ibd,time)

  implicit none
#include "SolidMechanics.h"
  logical, intent(in) :: restart
  integer, intent(in) :: ibd
  real, intent(in)    :: time

  !Local vars:
  character(len=6) :: str_ibd

  write(str_ibd,"(I6.6)") ibd

  ! When WRITEFORCES is defined --> for rigids , if clean start : we rewrite force and momt
  ! export files.
#ifdef WRITEFORCES
  write(*,*) 'INITIALIZE FORCES FILES, bod=',ibd
  if (.not. restart) then
     open(unit=113,file='./IOData/force.'//str_ibd//'.res',form='formatted', &
          status='replace')
     close(113)
     open(unit=113,file='./IOData/momt.'//str_ibd//'.res',form='formatted', &
          status='replace')
     close(113)
  endif
#endif

  ! When WRITESTATES is defined --> for rigids, if clean start : rewrite positions, vel and accel 
  ! export files.
  ! Translational vars
#ifdef WRITESTATES
  write(*,*) 'INITIALIZE STATES FILES, bod=',ibd
  if (.not. restart) then
     open(unit=113,file='./IOData/posvelacc_x.'//str_ibd//'.res',form='formatted', &
          status='replace')
     close(113)
  endif
#endif
  ! Orientation vars
#ifdef WRITESTATES
  if (.not. restart) then
     open(unit=113,file='./IOData/posvelacc_ang.'//str_ibd//'.res',form='formatted', &
          status='replace')
     close(113)
  endif
#endif
  ! Angular velocities 
#ifdef WRITESTATES
  if (.not. restart) then
     open(unit=113,file='./IOData/posvelacc_w.'//str_ibd//'.res',form='formatted', &
          status='replace')
     close(113)
  endif
#endif

  return

end subroutine sm_ioInit_rigid
