!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/GenAlpha/sm_GenAlpha_readCheckpoint
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
!!  Read in the sm_chkpt_ga.$(ibd).$(#).h5 file and set the values of
!!        qn, qdn, qddn, Hsn, Qsn
!!
!!***
#include "SolidMechanics.h"

subroutine sm_GenAlpha_readCheckpoint(ibd)
  use SolidMechanics_data, only: sm_structure, sm_BodyInfo
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use HDF5
  implicit none

  ! IO Varialbes
  integer, intent(in) :: ibd

  ! Internal Variables
  type(sm_structure),pointer :: body
  integer :: neq, ndofs, sm_checkpointFileNumber, h5err
  INTEGER(HID_T) :: file, dset, space
  INTEGER(HSIZE_T), DIMENSION(3) :: dims, chk_dims, chk_maxdims ! size read/write buffer
  logical :: link_exists
  character(LEN=100) :: filename

  ! set local pointer to body
  body => sm_bodyinfo(ibd)

  ! Set local value of the checkpoint number 
  call RuntimeParameters_get("checkPointFileNumber", sm_checkpointFileNumber)
  
  write(*,*) ''
  write(*,'(A,I4,A)',advance='no') 'Load GenAlpha CheckPt for body=',ibd,'...'
    
  !
  ! Make filename
  !
  write(filename,"(A,I1,A,I4.4,A)") 'sm_chkpt_ga.',ibd,'.',sm_checkpointFileNumber,'.h5'
  write(*,'(A)') trim(filename)

  !
  ! Check if File exists
  !
  inquire( file=trim(filename), exist= link_exists )
  if( .not. link_exists ) then
     call Driver_abortFlash("SM GenAlpha checkpoint file does not exist.")
  endif

  !
  ! Open the file
  !
  CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file, h5err) !H5F_ACC_RDONLY_F read only, _F fortran version
  call sm_io_checkHdfErr(h5err, 'failure to open body file')

  ! Set dims
  ndofs = body%ndofs
  neq   = body%neq

  !--------------------------------------------------------------
  ! qn(1:ndofs)
  !
  dims = (/ ndofs, 1, 1 /)

  ! Read in qn, it must be there
  call h5lexists_f( file, "/qn", link_exists, h5err )
  if( .not. link_exists ) then
     call Driver_abortFlash("qn not found")
  end if

  call h5dopen_f(file, "/qn", dset, h5err )

  ! check to make sure that qn is sized ndofs
  call H5Dget_space_f(dset, space, h5err)
  call H5Sget_simple_extent_dims_f(space, chk_dims, chk_maxdims, h5err) 
  call H5Sclose_f(space, h5err)
      
  if( int(chk_dims(1)) /= ndofs ) then
     call Driver_abortFlash("CheckPt size in file is not correct")
  end if

  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qn, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/qn')
  call h5dclose_f(dset, h5err)

  !--------------------------------------------------------------
  ! qdn(1:ndofs)
  !
  dims = (/ ndofs, 1, 1 /)

  ! Read in qn, it must be there
  call h5lexists_f( file, "/qdn", link_exists, h5err )
  if( .not. link_exists ) then
     call Driver_abortFlash("qdn not found")
  end if

  call h5dopen_f(file, "/qdn", dset, h5err )

  ! check to make sure that qdn is sized ndofs
  call H5Dget_space_f(dset, space, h5err)
  call H5Sget_simple_extent_dims_f(space, chk_dims, chk_maxdims, h5err) 
  call H5Sclose_f(space, h5err)
      
  if( int(chk_dims(1)) /= ndofs ) then
     call Driver_abortFlash("CheckPt size in file is not correct")
  end if

  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qdn, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/qdn')
  call h5dclose_f(dset, h5err)

  !--------------------------------------------------------------
  ! qddn(1:ndofs)
  !
  dims = (/ ndofs, 1, 1 /)

  ! Read in qddn, it must be there
  call h5lexists_f( file, "/qddn", link_exists, h5err )
  if( .not. link_exists ) then
     call Driver_abortFlash("qddn not found")
  end if

  call h5dopen_f(file, "/qddn", dset, h5err )

  ! check to make sure that qddn is sized ndofs
  call H5Dget_space_f(dset, space, h5err)
  call H5Sget_simple_extent_dims_f(space, chk_dims, chk_maxdims, h5err) 
  call H5Sclose_f(space, h5err)
      
  if( int(chk_dims(1)) /= ndofs ) then
     call Driver_abortFlash("CheckPt size in file is not correct")
  end if

  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%qddn, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/qddn')
  call h5dclose_f(dset, h5err)


  if( neq > 0 ) then

     !--------------------------------------------------------------
     ! Hsn(1:neq)
     !
     dims = (/ neq, 1, 1 /)

     ! Read in qddn, it must be there
     call h5lexists_f( file, "/Hs", link_exists, h5err )
     if( .not. link_exists ) then
        call Driver_abortFlash("Hs not found")
     end if

     call h5dopen_f(file, "/Hs", dset, h5err )

     ! check to make sure that Hs is size neq
     call H5Dget_space_f(dset, space, h5err)
     call H5Sget_simple_extent_dims_f(space, chk_dims, chk_maxdims, h5err) 
     call H5Sclose_f(space, h5err)

     if( int(chk_dims(1)) /= neq ) then
        call Driver_abortFlash("CheckPt size in file is not correct")
     end if

     call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Hsn, dims, h5err)
     call sm_io_checkHdfErr(h5err, '/Hs')
     call h5dclose_f(dset, h5err)

     !--------------------------------------------------------------
     ! Qsn(1:neq)
     !
     dims = (/ neq, 1, 1 /)

     ! Read in qddn, it must be there
     call h5lexists_f( file, "/Qs", link_exists, h5err )
     if( .not. link_exists ) then
        call Driver_abortFlash("Qs not found")
     end if

     call h5dopen_f(file, "/Qs", dset, h5err )

     ! check to make sure that Qs is size neq
     call H5Dget_space_f(dset, space, h5err)
     call H5Sget_simple_extent_dims_f(space, chk_dims, chk_maxdims, h5err) 
     call H5Sclose_f(space, h5err)

     if( int(chk_dims(1)) /= neq ) then
        call Driver_abortFlash("CheckPt size in file is not correct")
     end if

     call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Qsn, dims, h5err)
     call sm_io_checkHdfErr(h5err, '/Qs')
     call h5dclose_f(dset, h5err)  

  end if

  !--------------------------------------------------------------
  ! Fluid Forces: Hs_pres(1:ndofs) and Hs_visc(1:ndofs)
  !
  dims = (/ ndofs, 1, 1 /)

  ! Read in Hs_pres
  call h5lexists_f( file, "/Hs_pres", link_exists, h5err )
  if( .not. link_exists ) call Driver_abortFlash("Hs_pres not found")
  call h5dopen_f(file, "/Hs_pres", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Hs_pres, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/Hs_pres')
  call h5dclose_f(dset, h5err)

  ! Read in Hs_visc
  call h5lexists_f( file, "/Hs_visc", link_exists, h5err )
  if( .not. link_exists ) call Driver_abortFlash("Hs_visc not found")
  call h5dopen_f(file, "/Hs_visc", dset, h5err )
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%Hs_visc, dims, h5err)
  call sm_io_checkHdfErr(h5err, '/Hs_visc')
  call h5dclose_f(dset, h5err)
  
  !------------------------------------------------------------------------
  ! Close the file
  !
  CALL h5fclose_f(file , h5err)
  call sm_io_checkHdfErr(h5err, 'failure to close checkpoint file')

  ! Copy the values of qn -> qi
  body%qi(1:ndofs)   = body%qn(1:ndofs)
  body%qdi(1:ndofs)  = body%qdn(1:ndofs)
  body%qddi(1:ndofs) = body%qddn(1:ndofs)

  ! Copy old values of forces
  body%Hs(1:neq) = body%Hsn(1:neq)
  body%Qs(1:neq) = body%Qsn(1:neq)
  
  write(*,'(A)') '    complete.'

  return

end subroutine sm_GenAlpha_readCheckpoint
