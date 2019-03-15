!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_potential
!!
!! NAME 
!!
!!     Gravity_potential
!!
!! SYNOPSIS
!!
!!  call Gravity_potential(integer(IN) :: blockCount,
!!                                     integer(IN) :: blockList(blockCount),
!!                            optional,integer(IN) :: potentialIndex)
!!
!! DESCRIPTION
!!
!!      This routine computes the gravitational potential on all
!!      blocks specified in the list, for the gravity implementations
!!      (i.e., various Poisson implementations), which make use of it
!!      in computing the gravitational acceleration.
!!
!!      Supported boundary conditions are isolated (0) and
!!      periodic (1).  The same boundary conditions are applied
!!      in all directions.  For some implementation of Gravity,
!!      in particular with Barnes-Hut tee solver, additional combinations
!!      of boundary conditions may be supported.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!   potentialIndex : If present, determines which variable in UNK to use
!!                    for storing the updated potential.  If not present,
!!                    GPOT_VAR is assumed.
!!                    Presence or absense of this optional dummy argument
!!                    also determines whether some side effects are enabled
!!                    or disabled; see discussion of two modes under NOTES
!!                    below.
!!
!! NOTES
!!
!!  Gravity_potential can operate in one of two modes:
!!  * automatic mode  - when called without the optional potentialIndex.
!!    Such a call will usually be made once per time step, usually
!!    from the main time advancement loop in Driver_evolveFlash.
!!    Various side effects are enabled in this mode, see SIDE EFFECT below.
!!
!!  * explicit mode  - when called with the optional potentialIndex.
!!    The potential is stored in the variable explicitly given, and
!!    side effects like saving the previous potential in GPOL_VAR
!!    and updating some sink particle state and properties are
!!    suppressed.
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  gravitational potential.  Invokes a solver (of the Poisson equation)
!!  if necessary. On return, if potentialIndex is not present,
!!     GPOT_VAR:  contains potential for the current simulation time.
!!     GPOL_VAR (if defined): contains potential at the previous simulation time.
!!  On return, if potentialIndex is present, the UNK variable given by
!!  potentialIndex contains the newly computed potential.
!!
!!  May affect other variables related to particle properties if particles
!!  are included in the simulation.  In particular,
!!     PDEN_VAR (if defined): may get updated to the current density from
!!                particles if particles have mass.
!!
!!  There are additional side effects if sink particles are used.
!!  These effects happen by the call to Particles_sinkAccelGasOnSinksAndSinksOnGas,
!!  which may update sink particle properties and additional UNK variables that store
!!  accelerations. The calls are only made in automatic mode.
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Gravity implementation.
!!  The following information is subject to change without notice.
!!  For the Multigrid implementation:
!!     ISLS_VAR (residual)
!!     ICOR_VAR (correction)
!!     IMGM_VAR (image mass)
!!     IMGP_VAR (image potential)
!!  For the Multipole implementation:
!!     (none)
!!
!!***

!!REORDER(4): solnVec

subroutine Gravity_potential( potentialIndex)


  use Gravity_data, ONLY : grav_poisfact, grav_temporal_extrp, grav_boundary, &
       grav_unjunkPden, &
       useGravity, updateGravity, grv_meshComm
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY: Particles_updateGridVar, &
       Particles_sinkAccelGasOnSinksAndSinksOnGas
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
       Grid_getTileIterator, Grid_releaseTileIterator, &
       Grid_notifySolnDataUpdate, &
       Grid_solvePoisson
  use Grid_tile,     ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer, intent(IN), optional :: potentialIndex


  real, POINTER, dimension(:,:,:,:) :: solnVec

  integer       :: ierr

  real          :: redshift=0, oldRedshift=0
  real          :: scaleFactor, oldScaleFactor
  real          :: invscale, rescale
  integer       :: lb
  integer       :: bcTypes(6)
  real          :: bcValues(2,6) = 0.
  integer       :: density
  integer       :: newPotVar
  logical       :: saveLastPot
  type(Grid_tile_t)     :: tileDesc
  type(Grid_iterator_t) :: itor

  integer :: i, j, k

  nullify(solnVec)

  saveLastPot = (.NOT. present(potentialIndex))
  if (present(potentialIndex)) then
     newPotVar = potentialIndex
  else
     newPotVar = GPOT_VAR
  end if

  lb=1

!!$  call Cosmology_getRedshift(redshift)
!!$  call Cosmology_getOldRedshift(oldRedshift)
!!$  
  scaleFactor = 1./(1.+redshift)
  oldScaleFactor = 1./(1.+oldRedshift)
  
  invscale = 1./scaleFactor**3

! Rescaling factor to try and keep initial guess at potential close to
! final solution (in cosmological simulations).  Source term in Poisson
! equation has 1/a(t)^3 in it; and in linear theory (Omega=1, matter dom.)
! the comoving density increases as a(t), so comoving peculiar potential
! (which is what we are calculating here) should vary as 1/a(t)^2.  For
! noncosmological simulations this has no effect, since oldscale = scale = 1.

  rescale = (oldScaleFactor/scaleFactor)**2

!=========================================================================

  if(.not.useGravity) return
  
  if(.not.updateGravity) return

  call Timers_start("gravity Barrier")
  call MPI_Barrier (grv_meshComm, ierr)
  call Timers_stop("gravity Barrier")

  call Timers_start("gravity")

  bcTypes = grav_boundary
  where (bcTypes == PERIODIC)
     bcTypes = GRID_PDE_BND_PERIODIC
  elsewhere (bcTypes == ISOLATED)
     bcTypes = GRID_PDE_BND_ISOLATED
  elsewhere (bcTypes == DIRICHLET)
     bcTypes = GRID_PDE_BND_DIRICHLET
  elsewhere (bcTypes == OUTFLOW)
     bcTypes = GRID_PDE_BND_NEUMANN
  end where
  bcValues = 0.
 
  if (grav_temporal_extrp) then
     
     call Driver_abortFlash("shouldn't be here right now")
     !call extrp_initial_guess( igpot, igpol, igpot )
     
  else
     ! DEV: Do not enable tiling unless you have figured out the correct
     !      index space for setting solnVec in all locations below.
     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnVec, CENTER)
#ifdef GPOL_VAR
        if (saveLastPot) then
           ! DEV: The previous implementation was computing the data on the GC as
           !      well.  However, limited testing showed that it worked with just
           !      setting the data on the interior.  We default to setting the GC
           !      for now.  TBD if setting on the interior is the actual intent.
           do       k = tileDesc%grownLimits(LOW, KAXIS), tileDesc%grownLimits(HIGH, KAXIS)
              do    j = tileDesc%grownLimits(LOW, JAXIS), tileDesc%grownLimits(HIGH, JAXIS)
                 do i = tileDesc%grownLimits(LOW, IAXIS), tileDesc%grownLimits(HIGH, IAXIS)
                    solnVec(GPOL_VAR,i,j,k) = solnVec(GPOT_VAR,i,j,k)
                 end do
              end do
           end do
        end if
#endif
        ! DEV: The previous implementation was computing the data on the GC as
        !      well.  However, limited testing showed that it worked with just
        !      setting the data on the interior.  We default to setting the GC
        !      for now.  TBD if setting on the interior is the actual intent.
        do       k = tileDesc%grownLimits(LOW, KAXIS), tileDesc%grownLimits(HIGH, KAXIS)
           do    j = tileDesc%grownLimits(LOW, JAXIS), tileDesc%grownLimits(HIGH, JAXIS)
              do i = tileDesc%grownLimits(LOW, IAXIS), tileDesc%grownLimits(HIGH, IAXIS)
                 solnVec(newPotVar,i,j,k) = solnVec(newPotVar,i,j,k) * rescale
              end do
           end do
        end do

        ! CTSS - We should also be storing the old sink particle accelerations:
#if defined(SGXO_VAR) && defined(SGYO_VAR) && defined(SGZO_VAR)
        if (saveLastPot) then   !... but only if we are saving the old potential - kW
           call Driver_abortFlash("[Gravity_potential] Not tested third!")
           ! If tiling is used here, we probably need to write this as an
           ! explicit loop nest over the tile's indices
           solnVec(SGXO_VAR,:,:,:) = solnVec(SGAX_VAR,:,:,:)
           solnVec(SGYO_VAR,:,:,:) = solnVec(SGAY_VAR,:,:,:)
           solnVec(SGZO_VAR,:,:,:) = solnVec(SGAZ_VAR,:,:,:)
        end if
#endif

        ! for direct acceleration calculation by tree solver, added by R. Wunsch
#if defined(GAOX_VAR) && defined(GAOY_VAR) && defined(GAOZ_VAR)
        if (saveLastPot) then 
           call Driver_abortFlash("[Gravity_potential] Not tested fourth!")
           ! If tiling is used here, we probably need to write this as an
           ! explicit loop nest over the tile's indices
           solnVec(GAOX_VAR,:,:,:) = solnVec(GACX_VAR,:,:,:)
           solnVec(GAOY_VAR,:,:,:) = solnVec(GACY_VAR,:,:,:)
           solnVec(GAOZ_VAR,:,:,:) = solnVec(GACZ_VAR,:,:,:)
        end if
#endif
        call tileDesc%releaseDataPtr(solnVec, CENTER)
        call itor%next()
     enddo
     call Grid_releaseTileIterator(itor)
#ifdef GPOL_VAR
     if (saveLastPot) call Grid_notifySolnDataUpdate( (/GPOL_VAR/) )
#endif

  endif

! Poisson is solved with the total density of PDEN_VAR + DENS_VAR 
  density=DENS_VAR
! This only gets called if there are active particles.
#ifdef PDEN_VAR
#ifdef MASS_PART_PROP
  call Particles_updateGridVar(MASS_PART_PROP, PDEN_VAR)
  if (.NOT. grav_unjunkPden) call Grid_notifySolnDataUpdate( (/PDEN_VAR/) )
  density = PDEN_VAR
#ifdef DENS_VAR

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnVec, CENTER)
     call Driver_abortFlash("[Gravity_potential] Not tested either!")
     ! If tiling is used here, we probably need to write this as an
     ! explicit loop nest over the tile's indices
     solnVec(density,:,:,:) = solnVec(density,:,:,:) + &
          solnVec(DENS_VAR,:,:,:)
     call tileDesc%releaseDataPtr(solnVec, CENTER)
     call itor%next()
  enddo
  call Grid_releaseTileIterator(itor)
#endif
#endif
#endif

  invscale=grav_poisfact*invscale
  call Grid_solvePoisson (newPotVar, density, bcTypes, bcValues, &
       invscale)
  call Grid_notifySolnDataUpdate( (/newPotVar/) )

! Un-junk PDEN if it exists and if requested.

#ifdef PDEN_VAR
  if (grav_unjunkPden) then
     density = PDEN_VAR
#ifdef DENS_VAR           
     call Driver_abortFlash("[Gravity_potential] Not tested I guess!")
     ! I (JO) added this line to create the iterator.  I don't know if it is
     ! correct or if this code has ever been called.
     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnVec, CENTER)
        ! If tiling is used here, we probably need to write this as an
        ! explicit loop nest over the tile's indices
        solnVec(density,:,:,:) = solnVec(density,:,:,:) + &
             solnVec(DENS_VAR,:,:,:)
        call tileDesc%releaseDataPtr(solnVec, CENTER)
        call itor%next()
     enddo
     call Grid_releaseTileIterator(itor)
     if (density .NE. PDEN_VAR) call Grid_notifySolnDataUpdate( (/density/) )
#endif
  end if
#endif

  if (.NOT. present(potentialIndex)) then
     ! Compute acceleration of the sink particles caused by gas and vice versa
     call Particles_sinkAccelGasOnSinksAndSinksOnGas()
  end if
  
  
#ifdef USEBARS
  call MPI_Barrier (grv_meshComm, ierr)
#endif  
  call Timers_stop ("gravity")
  
  return
end subroutine Gravity_potential
