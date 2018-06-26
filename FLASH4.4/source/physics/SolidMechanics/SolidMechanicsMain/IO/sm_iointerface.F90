Module sm_iointerface
  
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"
  
  interface 
     subroutine sm_ioReadSolid(ibd)    
       implicit none
       integer, intent(IN) :: ibd
     end subroutine sm_ioReadSolid
  end interface
  
  interface 
     subroutine sm_io_checkHdfErr(hdferr, info_message)    
       implicit none
       integer, intent(in) :: hdferr
       CHARACTER, intent(in) :: info_message*(*)
     end subroutine sm_io_checkHdfErr
  end interface
  
  interface 
     subroutine sm_ioWriteSnapshot(ibd,idx_Snapshot)
       implicit none
       integer, intent(in) :: ibd, idx_Snapshot
     end subroutine sm_ioWriteSnapshot
  end interface
  
  interface 
     subroutine sm_ioRead_3DFlexible(ibd,file)
       use hdf5
       implicit none
       integer, intent(in) :: ibd
       INTEGER(HID_T), intent(in) :: file
     end subroutine sm_ioRead_3DFlexible
  end interface
  
  interface
     subroutine sm_ioRead_rbc(ibd,file)
       use hdf5
       implicit none
       integer, intent(IN) :: ibd
       integer(HID_T), intent(in) :: file
     end subroutine sm_ioRead_rbc
  end interface

  interface
     subroutine sm_ioRead_rigid(ibd,file)
       use hdf5
       implicit none
       integer, intent(IN) :: ibd
       integer(HID_T), intent(in) :: file
     end subroutine sm_ioRead_rigid
  end interface

  interface
     subroutine sm_ioRead_3DFlexible_loadIC(ibd)    
       implicit none
       integer, intent(IN) :: ibd
     end subroutine sm_ioRead_3DFlexible_loadIC
  end interface
  
  interface         
     subroutine sm_ioWrite_rbc(flag,mype,time,dt,istep,count, &
          timer,blockList,blockCount,firstfileflag) 
       implicit none 
       
       ! Argument list
       integer, intent(IN) ::mype,istep,count,firstfileflag
       real, intent(IN) :: time,dt,timer
       integer, intent(IN), optional :: flag
       integer, intent(IN) :: blockList,blockCount
      end subroutine sm_ioWrite_rbc
  end interface
  
  
  interface
       subroutine sm_io_rbcinit()
       implicit none

       end subroutine sm_io_rbcinit
  end interface

   interface 
     subroutine sm_xdmfwrite(filebase,ibd,mype,time,dt,istep,count, &
                          timer,firstfileflag) 
       implicit none
  
       ! Argument list
       character(len=100),intent (IN) :: filebase
       integer, intent(IN) :: mype,istep,count,firstfileflag,ibd
       real, intent(IN)  :: time,timer,dt  

   	end subroutine sm_xdmfwrite
   end interface

   interface 
	 subroutine outtoParaview(filename,mype,time,dt,istep,count, &
              timer,blockList,blockCount,firstfileflag)
           implicit none
           integer, intent(in)   :: mype,istep,count,firstfileflag
           integer, intent(in)   :: blockCount
           integer, intent(in)   :: blockList(MAXBLOCKS)
           real, intent(in)      :: time,dt,timer
           character(len=100),intent(in)  :: filename
         end subroutine outtoParaview
   end interface
   
   interface
      subroutine sm_ioWriteParticles(ibd,idx_Snapshot)
        implicit none  
        integer, intent(IN) :: ibd, idx_Snapshot
      end subroutine sm_ioWriteParticles
   end interface
      
   interface
      subroutine sm_ioInit_3DFlexible(restart,ibd,time)
        implicit none
        logical, intent(in) :: restart
        integer, intent(in) :: ibd
        real, intent(in)    :: time
      end subroutine sm_ioInit_3DFlexible
   end interface

   interface
      subroutine sm_ioInit_rbc(restart,ibd,time)
        implicit none
        logical, intent(in) :: restart
        integer, intent(in) :: ibd
        real, intent(in)    :: time
      end subroutine sm_ioInit_rbc
   end interface

   interface
      subroutine sm_ioInit_rigid(restart,ibd,time)
        implicit none
        logical, intent(in) :: restart
        integer, intent(in) :: ibd
        real, intent(in)    :: time
      end subroutine sm_ioInit_rigid
   end interface

   interface
      subroutine sm_ioInit(restart)
        implicit none
        logical, intent(in) :: restart
      end subroutine sm_ioInit
   end interface

   interface
      subroutine sm_iouttotecplot(nstep,time,dt,count)
        implicit none
        integer, intent(in) :: nstep,count
        real, intent(in)    :: time, dt
      end subroutine sm_iouttotecplot
   end interface

   interface
      subroutine sm_ioWriteStates_rigid(ibd,nstep,time)
        implicit none
        integer, intent(in) :: ibd,nstep
        real, intent(in)    :: time
      end subroutine sm_ioWriteStates_rigid
   end interface

   interface
      subroutine sm_ioWriteStates_3DFlexible(ibd,nstep,time)
        implicit none
        integer, intent(in) :: ibd,nstep
        real, intent(in)    :: time
      end subroutine sm_ioWriteStates_3DFlexible
   end interface

   interface
      subroutine sm_ioWriteStates_rbc(ibd,nstep,time)
        implicit none
        integer, intent(in) :: ibd,nstep
        real, intent(in)    :: time
      end subroutine sm_ioWriteStates_rbc
   end interface

   interface
      subroutine sm_ioWriteStates()
      end subroutine sm_ioWriteStates
   end interface


end Module sm_iointerface
