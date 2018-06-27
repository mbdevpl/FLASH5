!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/sm_ioWriteParticles
!!
!! NAME
!!  
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!!  
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_ioWriteParticles(ibd,idx_Snapshot)
   
  use SolidMechanics_Data, only : sm_structure,sm_BodyInfo
  use gr_sbData,           only : gr_sbBodyInfo
  use sm_iointerface,      only : sm_io_checkHdfErr
  use Driver_Interface,    only : Driver_getSimTime, Driver_getNstep
  USE HDF5
        
  implicit none
        
  ! IO Variables
  integer, intent(IN) :: ibd, idx_Snapshot

  ! Internal Variables
  INTEGER :: h5err, h5_outfile_id,group_id
  INTEGER(HID_T) :: dset_id, space_id
  INTEGER(HSIZE_T), DIMENSION(2) :: dims ! size read/write buffer
  character(len=100) ::   out_file_name
  real :: simulationTime
  integer :: totalPart, nel, nstep
  type(sm_structure), pointer :: body

  body => sm_BodyInfo(ibd)
  totalPart = gr_sbBodyInfo(ibd)%totalPart
  nel = sm_bodyInfo(ibd)%ws_nel

  write(*,'(A,I4,A)',advance='no') 'Writing particles for body=',ibd,'...'

  !
  ! Start writing a new hdf5 output file:
  !
  write(out_file_name,"(A,I5.5,A,I7.7,A)") 'IOData/sm_particles.',ibd,'.',idx_Snapshot,'.h5'
  call H5Fcreate_f(trim(out_file_name),H5F_ACC_TRUNC_F, h5_outfile_id, h5err)
  call sm_io_checkHDFerr(h5err, 'failure to open snapshots file')

  ! Change this if you want to write to a specific group
  group_id = h5_outfile_id

  !
  ! Write out the set of element nXi and nEta
  !
  dims = (/nel,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  ! nXi
  call H5Dcreate_f(group_id, 'nXi', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, sm_bodyInfo(ibd)%ws_nXi(1:nel), dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! nEta
  call H5Dcreate_f(group_id, 'nEta', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, sm_bodyInfo(ibd)%ws_nEta(1:nel), dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! close dataspace
  call H5Sclose_f(space_id, h5err)

  !
  ! Particle array size
  !
  dims = (/1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  ! Total number of particles
  call H5Dcreate_f(group_id, 'nParticles', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, totalPart, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! Total number of properties (redundant since it's in Flash.h)
  call H5Dcreate_f(group_id, 'NPART_PROPS', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, NPART_PROPS, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! close dataspace
  call H5Sclose_f(space_id, h5err)

  !
  ! Particle data
  ! 
  dims = (/ NPART_PROPS , totalPart /)
  call H5Screate_simple_f(2, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'particles', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  if( totalPart > 0 ) then
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,                               &
                     gr_sbBodyInfo(ibd)%particles(1:NPART_PROPS, 1:totalPart), &
                     dims, h5err)
  end if
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)

  !
  ! Write out the current time (for completeness if at different time from sm_results.*.h5)
  !
  call Driver_getSimTime( simulationTime )
  dims = (/1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'time', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, simulationTime, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)

  !
  ! Write out the current nstep
  !
  call Driver_getNstep( nstep )
  dims = (/1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'nstep', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, nstep, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)

  !
  ! Close the group
  !
  if( group_id /= h5_outfile_id ) then
     call H5Gclose_f (group_id, h5err)    
  end if

  !
  ! close the file
  !
  call H5Fclose_f(h5_outfile_id, h5err)    
  call sm_io_checkHDFerr(h5err,'failure to close snapshots file')

  ! complete.
  write(*,'(A)') 'complete'
  
end subroutine sm_ioWriteParticles
