!     
! File:   sm_pk_ReadHDF5.F90
! Author: tim
!
! 
subroutine sm_pk_ReadHDF5( infile_string )
      
      use sm_pk_data, only: sm_pk_info, sm_pk_NumKinematics
      use sm_iointerface, only: sm_io_checkHdfErr
      use SolidMechanics_data, only : sm_MeshMe
      use HDF5
      implicit none
#include "constants.h"      

      ! IO variables
      character, intent(in) :: infile_string*(*)
    
      ! internal variables   
      INTEGER :: h5err, i
      INTEGER(HID_T) :: file, dset
      INTEGER(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer
      character(len=100) :: read_str
    
      if(sm_meshMe .eq. MASTER_PE) write(*,'(A)',advance='no') 'Loading Kinematics Info...'
    
      !
      ! Open file
      !
      CALL h5fopen_f( trim( infile_string ), H5F_ACC_RDONLY_F, file, h5err)
      call sm_io_checkhdferr(h5err, 'failure to open kinematics file')
    
      ! read in number of kinematics sets to read
      dims = (/1,1,1/)
      CALL h5dopen_f (file, "/num_kinematics", dset, h5err)
      call sm_io_checkhdferr(h5err, '/num_kinematics')
      CALL h5dread_f(dset, H5T_NATIVE_INTEGER, sm_pk_NumKinematics, dims, h5err)
      CALL h5dclose_f(dset , h5err)

      ! Allocate the sm_pk_info
      allocate( sm_pk_info( sm_pk_NumKinematics ) )
      
      do i=1,sm_pk_NumKinematics
          ! Read in kinematics choice flag:
          dims = (/1,1,1/)
          write(read_str,"(A,I3.3,A)") '/kine',i,'/function_flag'
          CALL h5dopen_f (file, trim(read_str), dset, h5err)
          call sm_io_checkhdferr(h5err, trim(read_str) )
          CALL h5dread_f(dset, H5T_NATIVE_INTEGER, sm_pk_info(i)%flag, dims, h5err)
          CALL h5dclose_f(dset , h5err)

          ! read in number of parameters
          dims = (/1,1,1/)
          write(read_str,"(A,I3.3,A)") '/kine',i,'/num_parameters'
          CALL h5dopen_f (file, trim(read_str), dset, h5err)
          call sm_io_checkhdferr(h5err, trim(read_str) )
          CALL h5dread_f(dset, H5T_NATIVE_INTEGER, sm_pk_info(i)%NumParams, dims, h5err)
          CALL h5dclose_f(dset , h5err)

          ! resize the data storage
          allocate( sm_pk_info(i)%params( sm_pk_info(i)%NumParams ) )

          ! read in list of parameters
          dims = (/sm_pk_info(i)%NumParams,1,1/)
          write(read_str,"(A,I3.3,A)") '/kine',i,'/parameters'
          CALL h5dopen_f (file, trim(read_str), dset, h5err)
          call sm_io_checkhdferr(h5err, trim(read_str) )
          CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, sm_pk_info(i)%params, dims, h5err)
          CALL h5dclose_f(dset , h5err)

      end do
    
      !
      ! Close the file:
      !
      CALL h5fclose_f(file , h5err)
      call sm_io_checkhdferr(h5err, 'failure to close kinematics file')
      if(sm_meshMe .eq. MASTER_PE) write(*,'(A)') 'complete.'

end subroutine sm_pk_ReadHDF5
