!!****if* source/IO/IOMain/hdf5/parallel/io_closeFile
!!
!! NAME
!!  io_closeFile
!!
!! SYNOPSIS
!!
!!  io_closeFile(integer(in) :: fileID)
!!
!! DESCRIPTION
!!  
!!  closes the hdf5 file      
!!
!! ARGUMENTS
!!
!!  fileID - integer file identifier
!!
!!***

subroutine io_closeFile( fileID)

  use io_intfTypesModule, ONLY : io_fileID_t
  implicit none

  integer(io_fileID_t), intent(in) :: fileID
  call io_h5close_file(fileID)

end subroutine io_closeFile
