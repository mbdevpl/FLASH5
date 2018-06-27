subroutine  xmlwrite( filename,mype,count,dr_nstep)
  use Driver_data, Only : dr_nend
  implicit none
  !integer, intent(IN) :: nstring
  character (len=100),intent(IN):: filename
  integer, intent(IN)   :: mype,count,dr_nstep
  
  character (len=100) :: filenamexml,filenameh5
  integer :: XMLUNIT=50
  
  filenamexml='sm_rbc_animate.xml';
  write(*,*) 'The max no. of steps',dr_nend
  write(*,*)'This is step no. ',dr_nstep, count, mype

  write(filenameh5,'(A14,I4.4,".",I4.4,".xmf")') TRIM(filename),count,mype 
  
  if (dr_nstep.le.1) then
      
     open(XMLUNIT,file=TRIM(filenamexml),status='UNKNOWN');
     
     write(XMLUNIT,'(a)') '<?xml version="1.0" ?>'
     write(XMLUNIT,'(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
     write(XMLUNIT,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
     write(XMLUNIT,'(a)') '<Domain>'
     write(XMLUNIT,'(a)') '<Grid GridType="Collection" CollectionType="Temporal">'
     write(XMLUNIT,'(a,a,a)') '<xi:include href="',TRIM(filenameh5),'" xpointer="xpointer(//Xdmf/Domain/Grid)" />'
     close(XMLUNIT)
     write(*,*) 'First time: and filename is ',filenameh5
  elseif  ((dr_nstep.gt.1).and.(dr_nstep.lt.dr_nend)) then
     
     open(XMLUNIT,file=TRIM(filenamexml),status='UNKNOWN',POSITION='APPEND');
     write(XMLUNIT,'(a,a,a)') '<xi:include href="',TRIM(filenameh5),'" xpointer="xpointer(//Xdmf/Domain/Grid)" />'
     close(XMLUNIT)
     write(*,*) 'Second time: and filename is ',filenameh5
  elseif (dr_nstep.eq.dr_nend) then
     open(XMLUNIT,file=TRIM(filenamexml),status='OLD',POSITION='APPEND');
     write(XMLUNIT,'(a,a,a)') '<xi:include href="',TRIM(filenameh5),'" xpointer="xpointer(//Xdmf/Domain/Grid)" />'
     write(XMLUNIT,'(a)')'</Grid>'
     write(XMLUNIT,'(a)')'</Domain>'
     write(XMLUNIT,'(a)')'</Xdmf>'
     close(XMLUNIT)
     write(*,*) 'last time: and filename is ',filenameh5
  end if
end subroutine xmlwrite
