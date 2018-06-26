!     
! File:   sm_ioRead_3DFlexible_loadIC.F90
! Author: tim
! Notes: dof storage in body%qn, qdn has already been allocated

subroutine sm_ioRead_3DFlexible_loadIC(ibd)

#include "SolidMechanics.h"
    
  use SolidMechanics_Data, only: sm_structure,sm_BodyInfo
  use sm_iointerface, only: sm_io_checkHdfErr
  use Driver_interface, only: Driver_abortFlash

  USE HDF5
    
  implicit none
    
  integer, intent(IN) :: ibd

  INTEGER :: h5err, ndofs
  INTEGER(HID_T) :: file, dset, space
  INTEGER(HSIZE_T), DIMENSION(3) :: dims, chk_dims, chk_maxdims ! size read/write buffer
  logical :: link_exists
  character(LEN=100) :: filename

  type(sm_structure), pointer :: body
    

  body => sm_BodyInfo(ibd)

  write(*,*) ''
  write(*,'(A,I4,A)',advance='no') '   Load ICs for body=',ibd,'...'
    
  !
  ! Make filename
  !
  write(filename,"(A,I1,A)") 'sm_body_IC.',ibd,'.h5'
  write(*,'(A)') trim(filename)

  !
  ! Check if File exists
  !
  inquire( file=trim(filename), exist= link_exists )
  if( .not. link_exists ) then
     call Driver_abortFlash("SM IC file does not exist.")
  endif

  !
  ! Open the file
  !
  CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file, h5err) !H5F_ACC_RDONLY_F read only, _F fortran version
  call sm_io_checkHdfErr(h5err, 'failure to open body file')

  !
  ! Set dims
  !
  ndofs = body%ndofs
  dims = (/ ndofs, 1, 1 /)

  !
  ! Read in qn, it must be there
  !
  call h5lexists_f( file, "/qn", link_exists, h5err )
  if( .not. link_exists ) then
     call Driver_abortFlash("IC: qn not found")
  end if

  call h5dopen_f(file, "/qn", dset, h5err )

  ! check to make sure that qn is sized ndofs
  call H5Dget_space_f(dset, space, h5err)
  call H5Sget_simple_extent_dims_f(space, chk_dims, chk_maxdims, h5err) 
  call H5Sclose_f(space, h5err)
      
  if( chk_dims(1) /= ndofs ) then
     call Driver_abortFlash("IC size in file is not correct")
  end if

  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qn, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/qn')
  call h5dclose_f(dset, h5err)

  !
  ! Read in qdn, if it is there
  !
  call h5lexists_f( file, "/qdn", link_exists, h5err )
  if( link_exists ) then

      call h5dopen_f(file, "/qdn", dset, h5err )

      ! check to make sure that qn is sized ndofs
      call H5Dget_space_f(dset, space, h5err)
      call H5Sget_simple_extent_dims_f(space, chk_dims, chk_maxdims, h5err) 
      call H5Sclose_f(space, h5err)
      
      if( chk_dims(1) /= ndofs ) then
         call Driver_abortFlash("IC size in file is not correct")
      end if

      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qdn, dims, h5err)
      call sm_io_checkHdfErr(h5err, '/qdn')
      call h5dclose_f(dset, h5err)
  else
      write(*,'(A)') '      * qdn not found->assumming zeros.'
      body%qdn(1:ndofs) = 0.
  end if

  ! Copy over qn->qi
  body%qi(1:ndofs) = body%qn(1:ndofs)
  
  !
  ! Close the file
  !
  CALL h5fclose_f(file , h5err)
  call sm_io_checkHdfErr(h5err, 'failure to close IC file')
  
  write(*,'(A)') '    complete.'
  
end subroutine sm_ioRead_3DFlexible_loadIC


