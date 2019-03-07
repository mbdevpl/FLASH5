!!****h* source/physics/Gravity/Gravity_interface
!!
!! This is the header file for the gravity module that defines its
!! public interfaces.
!!***

Module Gravity_interface

#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

  interface Gravity_accel
     subroutine Gravity_accel (component,accelIndex, potentialIndex)
       integer, INTENT(in) ::  component
       integer, intent(in) :: accelIndex
       integer, intent(IN), optional :: potentialIndex
     end subroutine Gravity_accel
  end interface

  interface Gravity_accelOneRow
     subroutine Gravity_accelOneRow (pos,sweepDir,tileDesc, numCells, grav, Uin,&
           varIndex, extraAccelVars)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t),intent(IN) :: tileDesc
       integer, intent(IN) :: sweepDir,numCells
       integer, dimension(2),INTENT(in) ::pos
       real, dimension(numCells),INTENT(inout) :: grav
       real,    POINTER,    OPTIONAL :: Uin(:,:,:,:)
       integer, intent(IN), optional :: varIndex 
       integer, intent(IN),OPTIONAL      :: extraAccelVars(MDIM)
     end subroutine Gravity_accelOneRow
  end interface

  interface Gravity_accelOneBlock
     subroutine Gravity_accelOneBlock (tileDesc, ngcellcomp, gvec, potentialIndex)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t) :: tileDesc
      integer, intent(in)             :: ngcellcomp
       real, dimension(:,:,:,:),intent(out)  :: gvec
       integer, intent(in),optional    :: potentialIndex
     end subroutine Gravity_accelOneBlock
  end interface


  interface Gravity_computeDt
     subroutine Gravity_computeDt (Uin, dt_grav, dt_minloc)
       real,intent(OUT)       ::  dt_grav
       real,dimension(:,:,:,:) :: Uin
       integer, intent(INOUT) :: dt_minloc(5)
     end subroutine Gravity_computeDt
  end interface

  interface Gravity_finalize
     subroutine Gravity_finalize()
     end subroutine Gravity_finalize
  end interface

  interface Gravity_init
     subroutine Gravity_init()
       
     end subroutine Gravity_init
  end interface

  interface Gravity_potential
     subroutine Gravity_potential(potentialIndex)
       integer, intent(IN), optional :: potentialIndex
     end subroutine Gravity_potential
  end interface

  interface Gravity_unitTest
     subroutine Gravity_unitTest( fileUnit, perfect)
       implicit none
       integer, intent(in) :: fileUnit
       logical, intent(out) :: perfect
     end subroutine Gravity_unitTest
  end interface
end Module Gravity_interface
