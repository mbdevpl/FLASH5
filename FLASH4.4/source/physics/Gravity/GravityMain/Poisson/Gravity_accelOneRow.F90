!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           Grid_tile_t(IN)  :: tileDesc,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on a mesh block.  A single vector component
!!  of the acceleration is computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID   -     The local block metadata
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!  extraAccelVars      -  if specified,  Variables from which extra accelerations
!!                         are taken. Used to identify the UNK variables
!!                         that contain sink-on-gas accelerations when
!!                         sink particles are used.
!!
!! NOTES
!!
!!  If certain variables declared by the sink particles inplementation are
!!  declared, it is assumed that sink particles are in use.
!!  The sets of variables to make this determination are
!!    o  those given by extraAccelVars   extraAccelVars if present;
!!    o  {SGXO_VAR, SGYO_VAR, SGZO_VAR}  if potentialIndex is GPOL_VAR;
!!    o  {SGAX_VAR, SGAY_VAR, SGAZ_VAR}  otherwise.
!!  If it is assumed that sink particles are in use, then the acceleration
!!  returned in the grav array will have the appropriate sink particle
!!  acceleration component added to the acceleration computed by differencing
!!  the potential variable given by potentialIndex.
!!***

!!REORDER(4): solnVec

#include "constants.h"


subroutine Gravity_accelOneRow(pos, sweepDir, tileDesc, lo, hi, grav, Uin, &
                               potentialIndex, extraAccelVars)
  use Gravity_data, ONLY : grv_defaultGpotVar
  use Grid_tile,    ONLY : Grid_tile_t

  implicit none

#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling per-block drift logging for this routine because it is called too much
#undef Grid_releaseBlkPtr
#endif

  integer,           intent(IN)                      :: pos(2)
  integer,           intent(IN)                      :: sweepDir
  type(Grid_tile_t), intent(IN)                      :: tileDesc
  integer,           intent(IN)                      :: lo
  integer,           intent(IN)                      :: hi
  real,              intent(INOUT)                   :: grav(lo:hi)
  real,                            POINTER, OPTIONAL :: Uin(:,:,:,:)
  integer,           intent(IN),            OPTIONAL :: potentialIndex
  integer,           intent(IN),            OPTIONAL :: extraAccelVars(MDIM)

  real,   POINTER :: solnVec(:,:,:,:)
  integer         :: i
  real            :: gpot(lo:hi)
  real            :: delxinv
  integer         :: potVar
  integer         :: sink_ax_index, sink_ay_index, sink_az_index
  real            :: deltas(1:MDIM)

  !==================================================

  nullify(solnVec)
  
  call tileDesc%deltas(deltas)
  
  call tileDesc%getDataPtr(solnVec, CENTER)

!! IF a variable index is explicitly specified, assume that as the potential
!! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
     potVar=potentialIndex
  else if (grv_defaultGpotVar > 0) then
     potVar=grv_defaultGpotVar
  else
     potVar=GPOT_VAR
  end if

  sink_ax_index = 0
  sink_ay_index = 0
  sink_az_index = 0

  if (present(extraAccelVars)) then
  ! select the current or the old sink acceleration containers
     sink_ax_index = extraAccelVars(1)
     sink_ay_index = extraAccelVars(2)
     sink_az_index = extraAccelVars(3)
  else if (potVar .eq. GPOT_VAR) then
#if defined(SGAX_VAR) && defined(SGAY_VAR) && defined(SGAZ_VAR)
     sink_ax_index = SGAX_VAR
     sink_ay_index = SGAY_VAR
     sink_az_index = SGAZ_VAR
#endif
  else if (potVar .eq. GPOL_VAR) then
#if defined(SGXO_VAR) && defined(SGYO_VAR) && defined(SGZO_VAR)
     sink_ax_index = SGXO_VAR
     sink_ay_index = SGYO_VAR
     sink_az_index = SGZO_VAR
#endif
  endif
!!$  if ((potVar .ne. GPOT_VAR) .and. (potVar .ne. GPOL_VAR)) then
!!$     print *, "Gravity_accelOneRow called with neither GPOT_VAR nor GPOL_VAR."
!!$     call Driver_abortFlash("Gravity_accelOneRow called with neither GPOT_VAR nor GPOL_VAR.")
!!$  endif

  grav(lo:hi) = 0.0

  !Get row of potential values and compute inverse of zone spacing  
  if (sweepDir == SWEEP_X) then                     ! x-direction

     delxinv = 1.0 / deltas(IAXIS)

     do i = lo, hi 
        gpot(i) = solnVec(potVar,i,pos(1),pos(2))
     end do

     ! acceleration due to sink particles
     if (sink_ax_index > 0) then
        do i = lo, hi
           grav(i) = solnVec(sink_ax_index,i,pos(1),pos(2))
        end do
     end if

  elseif (sweepDir == SWEEP_Y) then                 ! y-direction

     delxinv = 1.0 / deltas(JAXIS)

     do i = lo, hi 
        gpot(i) = solnVec(potVar,pos(1),i,pos(2))
     end do

     ! acceleration due to sink particles
     if (sink_ay_index > 0) then
        do i = lo, hi
           grav(i) = solnVec(sink_ay_index,pos(1),i,pos(2))
        end do
     end if

  else                                            ! z-direction
     
     delxinv = 1.0 / deltas(KAXIS)

     do i = lo, hi 
        gpot(i) = solnVec(potVar,pos(1),pos(2),i)
     end do

     ! acceleration due to sink particles
     if (sink_az_index > 0) then
        do i = lo, hi 
           grav(i) = solnVec(sink_az_index,pos(1),pos(2),i)
        end do
     end if

  endif
  
  !-------------------------------------------------------------------------------
  
  !               Compute gravitational acceleration
  
  
  !**************** first-order differences
  !                 preserves conservation
  
  delxinv = 0.5e0 * delxinv
  
  do i = lo+1, hi-1
     grav(i) = grav(i) + delxinv * (gpot(i-1) - gpot(i+1))
  enddo
  
  grav(lo) = grav(lo+1)     ! this is invalid data - must not be used
  grav(hi) = grav(hi-1)
  
  call tileDesc%releaseDataPtr(solnVec, CENTER)
end subroutine Gravity_accelOneRow

