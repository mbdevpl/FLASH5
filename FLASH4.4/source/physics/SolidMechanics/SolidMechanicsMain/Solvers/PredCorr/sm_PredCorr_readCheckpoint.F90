!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/PredCorr/sm_PredCorr_readCheckpoint
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
!!  Read in the sm_chkpt_pc.$(ibd).$(#).h5 file and set the values of    
!!       pcmethod
!!       pciter
!!       body%qn(1:ndofs)
!!       body%qdn(1:ndofs)
!!       body%  qms(1:neq, -integ%pcmethod:-1)
!!       body% qdms(1:neq, -integ%pcmethod:-1)
!!       body%qddms(1:neq, -integ%pcmethod:-1)
!!       integ%vardt(-integ%pcmethod:0) )
!!
!!***
#include "SolidMechanics.h"

subroutine sm_PredCorr_readCheckpoint(ibd)
  use SolidMechanics_data, only: sm_structure, sm_BodyInfo
  use sm_PredCorr_data,    only: sm_PredCorr_type, sm_PredCorr_Info
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use HDF5
  implicit none

  ! IO Varialbes
  integer, intent(in) :: ibd

  ! Internal Variables
  type(sm_structure),     pointer   :: body
  type(sm_PredCorr_type), pointer   :: integ
  real, allocatable, dimension(:)   :: buffer1
  real, allocatable, dimension(:,:) :: buffer2
  integer :: neq, ndofs, sm_checkpointFileNumber, pcmethod_old, pciter_old, pcmethod
  integer :: h5err
  INTEGER(HID_T) :: file, dset
  INTEGER(HSIZE_T), DIMENSION(3) :: dims, chk_dims, chk_maxdims ! size read/write buffer
  logical :: link_exists
  character(LEN=100) :: filename

  ! set local pointers
  body  => sm_bodyinfo(ibd)
  integ => sm_PredCorr_Info(ibd)
  ndofs = body%ndofs
  neq   = body%neq

  ! Set local value of the checkpoint number 
  call RuntimeParameters_get("checkPointFileNumber", sm_checkpointFileNumber)  

  !------------------------------------------------------------------------
  ! Open the file
  !
  write(*,'(A,I4,A)',advance='no') 'Load PredCorr CheckPt for body=',ibd,'...'
    
  ! Make filename
  write(filename,"(A,I5.5,A,I4.4,A)") 'sm_chkpt_pc.',ibd,'.',sm_checkpointFileNumber,'.h5'
  write(*,'(A)') trim(filename)

  ! Check if File exists
  inquire( file=trim(filename), exist= link_exists )
  if( .not. link_exists ) then
     call Driver_abortFlash("SM PredCorr checkpoint file does not exist.")
  endif

  ! Open the file
  CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file, h5err)
  call sm_io_checkHdfErr(h5err, 'failure to open body file')

  !------------------------------------------------------------------------
  ! Compare pcmethod and pciter
  !
  dims = (/1,1,1/)
  
  ! pcmethod
  CALL h5dopen_f (file, "/pcmethod", dset, h5err) 
  call sm_io_checkHdfErr(h5err, '/pcmethod') 
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, pcmethod_old, dims, h5err)
  CALL h5dclose_f(dset , h5err)
 
  ! pciter
  CALL h5dopen_f (file, "/pciter", dset, h5err) 
  call sm_io_checkHdfErr(h5err, '/pciter') 
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, pciter_old, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  pcmethod = min( pcmethod_old, integ%pcmethod )
  integ%pciter = min( integ%pcmethod, pciter_old )

  !--------------------------------------------------------------
  ! qn, qdn, qddn
  !
  dims = (/ ndofs, 1, 1 /)

  ! qn:
  call h5dopen_f(file, "/qn", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qn, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/qn')
  call h5dclose_f(dset, h5err)

  ! qdn:
  call h5dopen_f(file, "/qdn", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qdn, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/qdn')
  call h5dclose_f(dset, h5err)

  ! qddn:
  call h5dopen_f(file, "/qddn", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qddn, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/qddn')
  call h5dclose_f(dset, h5err)

  !--------------------------------------------------------------
  ! qms, qdms, qddms
  !
  if (neq .gt. 0) then

    allocate( buffer2(neq, -pcmethod_old:-1) )
    dims = (/ neq, pcmethod_old, 1/)
  
    ! qms
    call h5dopen_f(file, "/qms", dset, h5err )
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, buffer2, dims, h5err)
    call sm_io_checkHdfErr(h5err, '/qms')
    call h5dclose_f(dset, h5err)
    body%qms(:,:) = 0.
    body%qms(1:neq, -pcmethod:-1) = buffer2(1:neq,-pcmethod:-1)

    ! qdms
    call h5dopen_f(file, "/qdms", dset, h5err )
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, buffer2, dims, h5err)
    call sm_io_checkHdfErr(h5err, '/qdms')
    call h5dclose_f(dset, h5err)
    body%qdms(:,:) = 0.
    body%qdms(1:neq, -pcmethod:-1) = buffer2(1:neq,-pcmethod:-1)

    ! qddms
    call h5dopen_f(file, "/qddms", dset, h5err )
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, buffer2, dims, h5err)
    call sm_io_checkHdfErr(h5err, '/qddms')
    call h5dclose_f(dset, h5err)
    body%qddms(:,:) = 0.
    body%qddms(1:neq, -pcmethod:-1) = buffer2(1:neq,-pcmethod:-1)

    ! deallocate buffer
    deallocate(buffer2)

  endif

  !--------------------------------------------------------------
  ! vartdt(-integ%pcmethod:0)
  !
  allocate( buffer1(-pcmethod_old:0) )
  dims = (/ pcmethod_old+1, 1, 1/)
  call h5dopen_f(file, "/vardt", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, buffer1, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/vardt')
  call h5dclose_f(dset, h5err)
  integ%vardt(:) = 0.
  integ%vardt(-pcmethod:0) = buffer1(-pcmethod:0)

  ! deallocate buffer
  deallocate(buffer1)

  !--------------------------------------------------------------
  ! Forces
  !

  if (neq .gt. 0) then

    dims = (/ neq, 1, 1/)
  
    ! Hs
    call h5dopen_f(file, "/Hs", dset, h5err )
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Hs, dims, h5err)
    call sm_io_checkHdfErr(h5err, '/Hs')
    call h5dclose_f(dset, h5err)

    ! Qs
    call h5dopen_f(file, "/Qs", dset, h5err )
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Qs, dims, h5err)
    call sm_io_checkHdfErr(h5err, '/Qs')
    call h5dclose_f(dset, h5err)

  endif

  ! Hs_pres, Hs_visc
  !
  dims = (/ ndofs, 1, 1 /)

  ! Hs_pres:
  call h5dopen_f(file, "/Hs_pres", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Hs_pres, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/Hs_pres')
  call h5dclose_f(dset, h5err)

  ! Hs_visc:
  call h5dopen_f(file, "/Hs_visc", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Hs_visc, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/Hs_visc')
  call h5dclose_f(dset, h5err)

  !------------------------------------------------------------------------
  ! Close the file
  !
  CALL h5fclose_f(file , h5err)
  call sm_io_checkHdfErr(h5err, 'failure to close checkpoint file')

  write(*,'(A)') '    complete.'

  return

end subroutine sm_PredCorr_readCheckpoint
