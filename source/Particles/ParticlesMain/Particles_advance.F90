!!****if* source/Particles/ParticlesMain/Particles_advance
!!
!! NAME
!!
!!  Particles_advance
!!
!! SYNOPSIS
!!
!!  Particles_advance(real(in) :: dtOld,
!!                    real(in) :: dtNew)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  Calls passive and active versions
!!  
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z} and VEL{X,Y,Z} properties of particles in the particles structure.
!!  Sorts particles in the particles structure by calling Grid_sortParticles.
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!***

!===============================================================================

subroutine Particles_advance (dtOld,dtNew)
  
  use Particles_data, ONLY: particles, pt_numLocal, pt_maxPerProc, useParticles, & 
       pt_meshMe, pt_typeInfo,&
       pt_indexList, pt_indexCount

  use Driver_interface, ONLY : Driver_abortFlash

  use pt_interface, ONLY: pt_updateTypeDS, pt_advanceRK
  use Grid_interface, ONLY : Grid_moveParticles, Grid_fillGuardCells, &
                             Grid_mapMeshToParticles
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
#include "GridParticles.h"
  
  real, INTENT(in)  :: dtOld, dtNew

  integer       :: i!,nstep,kk
  integer       ::p_begin,p_end
  logical,parameter :: regrid=.false.
  logical,save      :: gcMaskLogged = .FALSE.
  integer       :: pfor,pbak, lostNow

!!------------------------------------------------------------------------------
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return

  ! Prepare guardcell data needed for particle interpolation.
  !
  ! Experimentation with passive particles (with the old way of advancing particles)
  ! has shown that at least 2 layers of guardcells need to be filled
  ! with updated data for vel[xyz] and density, in order to get the
  ! same results as for a full guardcell fill, when using native grid interpolation. - KW
  ! With "monotonic" interpolation, even more layers are needed. - KW
  p_begin=1
  p_end=pt_numLocal
#ifndef PRNT_PART_PROP
!!$  if (pt_reduceGcellFills) then
!!$     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,unitReadsMeshDataOnly=.true.)
!!$  else
     call Grid_fillGuardCells( CENTER, ALLDIR)
!!$          maskSize=pt_gcMaskSizeForAdvance,mask=pt_gcMaskForAdvance,&
!!$          doLogMask=.NOT.gcMaskLogged)
!!$  end if
#endif

!!$  els
  ! Now update the pt_typeInfo data structure
  !! ?? Why is call to pt_updateTypeDS required after Grid_sortParticles??
!!$  call pt_updateTypeDS(particlesPerBlk)
  

  !! Now do actual movement, advance particles in time.
  !! This implements the option of picking different
  !! integration methods for different particle types.
  i=1
  select case(pt_typeInfo(PART_ADVMETHOD,i))
  case(RUNGEKUTTA) 
     call pt_advanceRK(dtOld,dtNew, p_begin,p_end,i)
!DevNote :: Following two options to be implemented later
!  case(MIDPOINT)
!     call pt_advanceMidpoint(dtOld,dtNew,p_begin,p_end,i)
!  case(EULER_TRA)
!     call pt_advanceEuler_passive(dtOld,dtNew,p_begin,p_end,i)
  case default
     call Driver_abortFlash("Particles_advance: Not a valid advance method. Please use RUNGEKUTTA method!")
  end select
  

#ifdef DEBUG_PARTICLES
  print*,' ready to move Particles'
#endif

#ifdef DEBUG_VPARTICLES
  do kk=1,pt_numLocal
     write(*,*)'local particle ',kk,'In proc ',int(particles(PROC_PART_PROP,kk)),'on blk=',int(particles(BLK_PART_PROP,kk))
     write(*,*)'Xpos= ',particles(POSX_PART_PROP,kk),'Ypos=',particles(POSY_PART_PROP,kk)
  enddo
#endif 
  
  ! Put the particles in the appropriate blocks if they've moved off
  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
       pt_indexList, pt_indexCount, regrid) 

  ! Now update the pt_typeInfo data structure
  !! ?? Why is call to pt_updateTypeDS required after Grid_sortParticles??
!!$  call pt_updateTypeDS(particlesPerBlk)

  
  ! If predictive routines are used, they will need to sort and prepare for the
  !  next time step.  Since sorting is so expensive, we suffer code duplication
  !  and do it in the pt_preparePassive routines.
  ! Many algorithms use the stub routines.

  gcMaskLogged = .TRUE.
  
  return

  !!-----------------------------------------------------------------------
end subroutine Particles_advance


