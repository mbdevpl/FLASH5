!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/GenAlpha/sm_GenAlpha_writeCheckpoint
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
!!  write the sm_chkpt_ga.$(ibd).$(#).h5 file and set the values of
!!        qn, qdn, qddn, Hs, Qs, dt, nstep, time
!!
!!***
#include "SolidMechanics.h"

subroutine sm_GenAlpha_writeCheckpoint(ibd, checkpt_num)
  use SolidMechanics_data, only: sm_structure, sm_BodyInfo
  use Driver_interface,    only: Driver_getSimTime, Driver_getNstep
  use sm_surf_interface,   only: sm_surf_currentFluidForce
  use HDF5
  implicit none

  ! IO Varialbes
  integer, intent(in) :: ibd
  integer, intent(in) :: checkpt_num

  ! Internal Variables
  type(sm_structure),pointer :: body
  real :: dt, time
  integer :: neq, ndofs, nstep, h5err
  INTEGER(HID_T) :: file, dset, space
  INTEGER(HSIZE_T), DIMENSION(3) :: dims
  logical :: link_exists
  character(LEN=100) :: filename
  real, allocatable, dimension(:) :: Hs_pres, Hs_visc

  ! set local pointer to body
  body => sm_bodyinfo(ibd)

  write(*,'(A,I4,A)',advance='no') 'write GenAlpha CheckPt for body=',ibd,'...'
    
  !
  ! Make filename
  !
  write(filename,"(A,I1,A,I4.4,A)") 'sm_chkpt_ga.',ibd,'.', checkpt_num,'.h5'
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
  ! qn, qdn, qddn
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

  ! qddn:
  call H5Dcreate_f(file, 'qddn', H5T_NATIVE_DOUBLE, space, dset, h5err)
  call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%qddn, dims, h5err)
  call H5Dclose_f(dset, h5err)

  ! close the space
  call H5Sclose_f(space, h5err)


  !------------------------------------------------------------- 
  ! For body not fully restrained
  ! Hs, Qs
  !
  if( neq > 0 ) then 
  
     dims = (/neq, 1, 1/)
     call H5Screate_simple_f(1, dims, space, h5err)

     ! Qs
     call H5Dcreate_f(file, 'Qs', H5T_NATIVE_DOUBLE, space, dset, h5err)
     call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%Qs, dims, h5err)
     call H5Dclose_f(dset, h5err)

     ! Hs
     call H5Dcreate_f(file, 'Hs', H5T_NATIVE_DOUBLE, space, dset, h5err)
     call H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, body%Hs, dims, h5err)
     call H5Dclose_f(dset, h5err)

      ! close the space
     call H5Sclose_f(space, h5err)

  end if

  !------------------------------------------------------------- 
  ! Fluid forces, at end of step
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

  return

end subroutine sm_GenAlpha_writeCheckpoint
