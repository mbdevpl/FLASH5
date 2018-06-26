!     
! File:   init.F95
! Author: tim
!
#include "constants.h"
#include "SolidMechanics.h"
#include "Flash.h"

subroutine sm_ioReadSolid(ibd)
  use SolidMechanics_Data, only : sm_MeshMe,sm_structure,sm_BodyInfo
  use sm_iointerface, only: sm_io_checkHdfErr, sm_ioRead_3DFlexible, sm_ioRead_rigid, sm_ioRead_rbc
  USE HDF5
  implicit none
    
  integer, intent(IN) :: ibd
  character(LEN=100)  :: filename

  INTEGER :: h5err
  INTEGER(HID_T),save :: file, dset
  INTEGER(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer

  type(sm_structure), pointer :: body

  body => sm_BodyInfo(ibd)

  !
  ! Open file
  !
  if (sm_meshMe==body%BodyMaster) then


     if (body%MetaBody .gt. CONSTANT_ZERO) then
       ! File number takes metabody index
       write(filename,"(A,I5.5,A)") 'sm_body.',body%MetaBody,'.h5'
     else
       ! File number takes body index
       write(filename,"(A,I5.5,A)") 'sm_body.',ibd,'.h5'
     endif

     CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file, h5err) !H5F_ACC_RDONLY_F read only, _F fortran version
     call sm_io_checkHdfErr(h5err, 'failure to open body file')

     !
     ! Check BodyType
     !
     dims = (/1,1,1/)
     CALL h5dopen_f (file, "BodyType", dset, h5err) ! dset handle to datase required
     call sm_io_checkHdfErr(h5err, 'BodyType') ! error check
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%BodyType, dims, h5err) !use handle
     CALL h5dclose_f(dset , h5err)

     ! Load and allocate based on BodyType
     select case( body%BodyType )

     case( BODYTYPE_RIGID ) 
        if (ibd .eq. 1) write(*,'(A)',advance='no') '... Rigid ...'
        call sm_ioRead_Rigid(ibd,file)

     case( BODYTYPE_2DFLEXIBLE )
        if (ibd .eq. 1) write(*,'(A)',advance='no') '... 2D Flexible ...'
        call Driver_abortFlash('BodyType not yet implemented')

     case( BODYTYPE_3DFLEXIBLE )
        if (ibd .eq. 1) write(*,'(A)',advance='no') '... 3D Flexible ...'
        call sm_ioRead_3DFlexible(ibd,file)

     case( BODYTYPE_RBC )
        if (ibd .eq. 1) write(*,'(A)',advance='no') '... RBC ...'
        call sm_ioRead_RBC(ibd,file)

     case default
        call Driver_abortFlash("BodyType not yet implemented")

     end select

     !
     ! Close the file:
     !
     call h5fclose_f(file, h5err)
     call sm_io_checkHdfErr(h5err, 'failure to close mesh file')

  end if

end subroutine sm_ioReadSolid
