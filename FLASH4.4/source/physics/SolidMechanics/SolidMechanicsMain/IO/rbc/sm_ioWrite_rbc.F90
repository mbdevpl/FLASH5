
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!          This subroutine write the RBC positions in ASCII format 
!          as Finite element data at each output step
!          written by Hussein Ezzeldin July 2010
!          future enhancments write the data in Binary format
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
subroutine sm_ioWrite_rbc(flag,mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag) 

  
#include "Flash.h"
#include "SolidMechanics.h"

  use SolidMechanics_data, ONLY : sm_structure, sm_bodyInfo, sm_NumBodies
  use SolidMechanics_rbc_data, Only : rbcplotOutputInterval
  use sm_iointerface, Only : outtoParaview
  implicit none 
      
  ! Argument list
  integer, intent(IN) ::mype,istep,count,firstfileflag
  real, intent(IN) :: time,dt,timer
  integer, intent(IN), optional :: flag
  integer, intent(IN) :: blockList(MAXBLOCKS),blockCount
  
  integer :: i,j,filetype=1;
  integer :: lstring
  character(len=100) :: filename,RBCname
  character(len=100) :: filename1,Fname
  type(sm_structure),pointer :: body
  !real,pointer :: NCVALUES, CCVALUES
  
!     Create the file name
  if (present(flag)) filetype=flag;
  
  
  write(*,*) 'The count at top of sm_ioWrite_RBC is: ',count
  select case (filetype) 
  case (WRITEPOS)
      !body => sm_BodyInfo(ibd)	  
      !write(filename,'(A13,I9.9,A4)') './IOData/RBC_',istep,'.plt'
      !write(*,*) 'RBC positions  are written to file: ',filename
     
      !open(unit=11,file=filename,status='unknown')
     
     !
      !write(11,*) 'TITLE = "RBC"'
      !write(11,*) 'VARIABLES = "X", "Y", "Z"'
     
     !     Writing coordinates
     !write(RBCname,'("RBC",I9.9)') istep
     
     !write(11,'("ZONE T=",A," DATAPACKING=POINT, NODES=",   &
         !! & I0,", ELEMENTS=",I0,", ZONETYPE=FETRIANGLE")') TRIM(RBCname), &
          !body%nnp, body%nele
     
     !do i=1,body%nnp
      !  write(11,'(3F14.5)')body%qms(3*i-2:3*i,1)
     !enddo
     
     !write(11,*)' '
     
    ! do i=1, body%nele
    !    write(11,*)body%Tri(1:3,i)
    ! enddo
     
    ! close(11)
     
  case (WRITEVEL)
     
  case (WRITEFORC)
!c$$$c$$$     body => sm_BodyInfo(ibd)
!c$$$c$$$     
!c$$$c$$$     write(filename1,'(A16,I9.9,A4)') './FORCES_',istep,'.dat'
!c$$$c$$$     write(*,*) 'Nodal Forces are written to file: ',TRIM(filename1)
!c$$$c$$$     open(unit=12,file=filename1,status='unknown')
!c$$$c$$$     write(12,*) 'TITLE = "F_NODES"'
!c$$$c$$$     write(12,*) 'VARIABLES = "FX", "FY", "FZ"'
!c$$$c$$$     
!c$$$c$$$     
!c$$$c$$$     write(Fname,'("F_",I9.9)') istep
!c$$$c$$$     
!c$$$c$$$     write(12,'("ZONE T=",A," DATAPACKING=POINT, NODES=",  &
!c$$$c$$$          & I0,", ELEMENTS=",I0,", ZONETYPE=FETRIANGLE")') TRIM(Fname), &
!c$$$c$$$          body%nnp, body%nele
!c$$$c$$$     
!c$$$c$$$     
!c$$$c$$$     do i=1,body%nnp
!c$$$c$$$        write(12,'(3F14.5)')body%qddms(3*i-2:3*i,1)
!c$$$c$$$     enddo
!c$$$c$$$     
!c$$$c$$$     write(12,*)' '
!c$$$c$$$     
!c$$$c$$$     do i=1, body%nele  
!c$$$c$$$        write(12,*)body%Tri(1:3,i)
!    enddo
     
  case (WRITEHDF5)

     filename='sm_rbc_output_';
     write(*,*) 'Calling outtoParaview'
     call outtoParaview(filename,mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)

     !call sm_xdmfwrite(filename,ibd,mype,time,dt,istep,count, &
     !                     timer,firstfileflag)
     !filename1='sm_rbc_animate.xml';
     write(*,*) 'firstfileflag',firstfileflag
     call xmlwrite ( filename,mype,count,istep)
     
     write(filename,'(A13,I9.9,A4)') './IOData/ContantData_',istep,'.txt'
  end select
  
  
end subroutine sm_ioWrite_rbc

