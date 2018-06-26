!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/PredCorr/sm_PredCorr_writeCheckpoint
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
!!  Write out the sm_chkpt_pc.$(ibd).$(#).h5 file from the values of    
!!       pcmethod
!!       pciter
!!       body%qn(1:ndofs)
!!       body%qdn(1:ndofs)
!!       body%Hs(1:neq)
!!       body%  qms(1:neq, -integ%pcmethod:-1)
!!       body% qdms(1:neq, -integ%pcmethod:-1)
!!       body%qddms(1:neq, -integ%pcmethod:-1)
!!       integ%vardt(-integ%pcmethod:0) )
!!       time
!!       nstep
!!
!!***
#include "SolidMechanics.h"

subroutine sm_PredCorr_writeCheckpoint(ibd, checkpt_num)
  use SolidMechanics_data, only: sm_structure, sm_BodyInfo
  use sm_PredCorr_data,    only: sm_PredCorr_type, sm_PredCorr_Info
  use sm_surf_interface,   only : sm_surf_currentFluidForce
  use Driver_interface,    only: Driver_getSimTime, Driver_getNstep
  use HDF5
  implicit none

  ! IO Varialbes
  integer, intent(in) :: ibd
  integer, intent(in) :: checkpt_num

  ! Internal Variables
  type(sm_structure),     pointer :: body
  type(sm_PredCorr_type), pointer :: integ
  real :: dt, time
  integer :: neq, ndofs, nstep, h5err
  INTEGER(HID_T) :: file, dset, space
  INTEGER(HSIZE_T), DIMENSION(3) :: dims
  logical :: link_exists
  character(LEN=100) :: filename
  real, allocatable, dimension(:) :: Hs_pres, Hs_visc


  ! set local pointer to body
  body => sm_bodyinfo(ibd)

  ! set integ pointer:
  integ => sm_PredCorr_info(ibd)

  write(*,'(A,I4,A)',advance='no') 'write PredCorr CheckPt for body=',ibd,'...'
    
  !
  ! Make filename
  !
  write(filename,"(A,I5.5,A,I4.4,A)") 'sm_chkpt_pc.',ibd,'.',checkpt_num,'.h5'
  write(*,'(A)') trim(filename)

  !
  ! Open the file
  !
  CALL h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file, h5err)
  call sm_io_checkHdfErr(h5err, 'failure to open chkpt file for writing')

  ! Set dims
  ndofs = body%ndofs
  neq   = body%neq

  !------------------------------------------------------------- 
  ! degrees of freedom
  ! qn, qdn
  !
  dims = (/ndofs, 1, 1/)
  call H5Screate_simple_f(1, dims, space, h5err)

  ! qn
  call H5Dcreate_f(file, 'qn', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%qn, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! qdn
  call H5Dcreate_f(file, 'qdn', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%qdn, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! qddn
  call H5Dcreate_f(file, 'qddn', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%qddn, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! close the space
  call H5Sclose_f(space, h5err)

  !------------------------------------------------------------- 
  ! qms, qdms, qddms
  !

  if ( neq .gt. 0 ) then
    dims = (/neq, integ%pcmethod, 1/)
    call H5Screate_simple_f(1, dims, space, h5err)

    ! qms
    call H5Dcreate_f(file, 'qms', H5T_NATIVE_DOUBLE, space, dset, h5err)
    call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%qms, dims, h5err)
    call H5Dclose_f(dset, h5err)

    ! qdms
    call H5Dcreate_f(file, 'qdms', H5T_NATIVE_DOUBLE, space, dset, h5err)
    call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%qdms, dims, h5err)
    call H5Dclose_f(dset, h5err)

    ! qddms
    call H5Dcreate_f(file, 'qddms', H5T_NATIVE_DOUBLE, space, dset, h5err)
    call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%qddms, dims, h5err)
    call H5Dclose_f(dset, h5err)

    ! close the space
    call H5Sclose_f(space, h5err)

  endif

  !------------------------------------------------------------- 
  ! Forces
  !
  if( neq .gt. 0 ) then

    dims = (/ neq, 1, 1/)
    call H5Screate_simple_f(1, dims, space, h5err)

    ! Hs
    call H5Dcreate_f(file, 'Hs', H5T_NATIVE_DOUBLE, space, dset, h5err)
    call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%Hs, dims, h5err)
    call H5Dclose_f(dset, h5err)

    ! Qs
    call H5Dcreate_f(file, 'Qs', H5T_NATIVE_DOUBLE, space, dset, h5err)
    call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%Qs, dims, h5err)
    call H5Dclose_f(dset, h5err)

    ! close the space
    call H5Sclose_f(space, h5err)

  endif

  ! Hs_pres, Hs_visc
  !
  allocate( Hs_pres(body%ndofs), Hs_visc(body%ndofs) )
  call sm_surf_currentFluidForce(ibd, body%ndofs, Hs_pres, Hs_visc)

  dims = (/ndofs, 1, 1/)
  call H5Screate_simple_f(1, dims, space, h5err)

  ! Hs_pres
  call H5Dcreate_f(file, 'Hs_pres', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, Hs_pres, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! Hs_visc
  call H5Dcreate_f(file, 'Hs_visc', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, Hs_visc, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! close the space
  call H5Sclose_f(space, h5err)



  !------------------------------------------------------------- 
  ! vardt, pcmethod, pciter
  !
  dims = (/integ%pcmethod+1, 1, 1/)
  call H5Screate_simple_f(1, dims, space, h5err)

  ! vardt
  call H5Dcreate_f(file, 'vardt', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, integ%vardt, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! close the space
  call H5Sclose_f(space, h5err)

  dims = (/1, 1, 1/)
  call H5Screate_simple_f(1, dims, space, h5err)
  ! PC method:
  call H5Dcreate_f(file, 'pcmethod', H5T_NATIVE_INTEGER, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_INTEGER, integ%pcmethod, dims, h5err)
  call H5Dclose_f(dset, h5err)  

  ! PC iter:
  call H5Dcreate_f(file, 'pciter', H5T_NATIVE_INTEGER, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_INTEGER, integ%pciter, dims, h5err)
  call H5Dclose_f(dset, h5err)  

  ! close the space
  call H5Sclose_f(space, h5err)
   
  !------------------------------------------------------------- 
  !  Write out scalars
  !
  dims = (/1,1,1/)
  call H5Screate_simple_f(1, dims, space, h5err)

  ! Simulation Time
  call Driver_getSimTime( time )
  call H5Dcreate_f(file, 'time', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, time, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! Simulation Nstep
  call Driver_getNstep( nstep )
  call H5Dcreate_f(file, 'nstep', H5T_NATIVE_INTEGER, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_INTEGER, nstep, dims, h5err)
  call H5Dclose_f(dset, h5err)


  ! close the space
  call H5Sclose_f(space, h5err)

  !
  ! close the file
  !
  call H5Fclose_f(file, h5err)  
  call sm_io_checkHDFerr(h5err,'failure to close checkpoint file')

  deallocate(Hs_pres, Hs_visc)

  return

end subroutine sm_PredCorr_writeCheckpoint
