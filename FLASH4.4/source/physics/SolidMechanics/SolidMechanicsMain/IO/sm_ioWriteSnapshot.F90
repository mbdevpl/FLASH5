!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/sm_ioWriteSnapshot
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

subroutine sm_ioWriteSnapshot(ibd,idx_Snapshot)
   
  use SolidMechanics_Data,  only: sm_structure,sm_BodyInfo
  use sm_iointerface,       only: sm_io_checkHdfErr
  use sm_integdata,         only: sm_integ_subiter_old
  use sm_assemble_interface,only: sm_assemble_COM
  use sm_surf_interface,    only: sm_surf_assembleFluidForce_toPoint
  use Driver_Interface,     only: Driver_getSimTime, Driver_getNstep
  USE HDF5
        
  implicit none

  ! IO Variables
  integer, intent(IN) :: ibd, idx_Snapshot

  ! Internal Variables
  INTEGER :: h5err, h5_outfile_id,group_id
  INTEGER(HID_T) :: dset_id, space_id
  INTEGER(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer
  character(len=100) ::   out_file_name
  real :: simulationTime
  real, dimension(NDIM) :: force_pres, force_visc, moment_pres, moment_visc
  integer :: nstep
  type(sm_structure), pointer :: body
    
  body => sm_BodyInfo(ibd)
      
  write(*,'(A,I4,A)',advance='no') 'Writing snapshot for body=',ibd,'...'

  !
  ! Start writing a new hdf5 output file:
  !
  write(out_file_name,"(A,I5.5,A,I7.7,A)") 'IOData/sm_results.',ibd,'.',idx_Snapshot,'.h5'
  call H5Fcreate_f(trim(out_file_name),H5F_ACC_TRUNC_F, h5_outfile_id, h5err)
  call sm_io_checkHDFerr(h5err, 'failure to open snapshots file')
  
  ! Change this if you want to write to a specific group
  group_id = h5_outfile_id
      
  !
  ! Write qn:
  !
  dims = (/body%ndofs,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'qn', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%qn, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  
  ! qdn:
  call H5Dcreate_f(group_id, 'qdn', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%qdn, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  
  ! qddn:
  call H5Dcreate_f(group_id, 'qddn', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%qddn, dims, h5err)
  call H5Dclose_f(dset_id, h5err)    
  
  ! close the space:
  call H5Sclose_f(space_id, h5err)
 

  if (body%neq .gt. 0) then 
  !
  ! Write out the current internal loading information
  !
  dims = (/body%neq,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'Qs', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%Qs, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)
  
  !
  ! Write out the current body loads
  !
  dims = (/body%neq,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'Hs', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%Hs, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)
  endif

  !
  ! Write out the current body fluid loads
  !
  dims = (/body%ndofs,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  ! Pressure
  call H5Dcreate_f(group_id, 'Hs_pres', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%Hs_pres, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! Visc.
  call H5Dcreate_f(group_id, 'Hs_visc', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%Hs_visc, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! Close the dataspace
  call H5Sclose_f(space_id, h5err)
  
  !
  ! Write out the current time
  !
  call Driver_getSimTime( simulationTime )
  dims = (/1,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'time', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, simulationTime, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)

  !
  ! Write out the current nstep
  !
  call Driver_getNstep( nstep )
  dims = (/1,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'nstep', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, nstep, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)
  
  !
  ! Number of subsets
  !
  dims = (/1,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'SubIt', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, sm_integ_subiter_old, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)
  
  !
  ! Center of Mass
  !
  ! update body%COM
  call sm_assemble_COM(ibd)
  dims = (/NDIM,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id, 'CenterOfMass', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%COM, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)
  
  !
  ! Total Force and Moment of Fluid
  !
  call sm_surf_assembleFluidForce_toPoint(ibd, body%COM, &
                                          force_pres, force_visc, moment_pres, moment_visc)
  dims = (/NDIM,1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)  
  ! Pressure: force
  call H5Dcreate_f(group_id, 'FluidForce_pres', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, force_pres, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! Pressure moment
  call H5Dcreate_f(group_id, 'FluidMoment_pres', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, moment_pres, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! Visc. force
  call H5Dcreate_f(group_id, 'FluidForce_visc', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, force_visc, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! Visc. moment
  call H5Dcreate_f(group_id, 'FluidMoment_visc', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, moment_visc, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  ! close dataspace
  call H5Sclose_f(space_id, h5err)
  
  ! Close the group
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
  
end subroutine sm_ioWriteSnapshot
