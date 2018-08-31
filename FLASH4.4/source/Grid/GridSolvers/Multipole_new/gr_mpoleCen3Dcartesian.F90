!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCen3Dcartesian
!!
!! NAME
!!
!!  gr_mpoleCen3Dcartesian
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen3Dcartesian (integer, intent (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 3D cartesian
!!  geometry. It computes the center of expansion for the multipoles for 3D cartesian
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                            integral (r * rho * rho  dr)
!!              Cen (x,y,z) = ----------------------------
!!                              integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                     sum cells (cell center r * cell mass * cell rho)
!!       Cen (x,y,z) = ------------------------------------------------
!!                             sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (x,y,z) has been determined, it is placed on the
!!  the nearest cell corner. The following is computed here:
!!
!!                  1) multipole expansion center (placed on nearest cell corner)
!!                  2) total mass (aborts, if <= 0)
!!                  3) the 'atomic' inner zone length (and its inverse)
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen3Dcartesian (idensvar)

  use Grid_data,         ONLY : gr_meshMe,  &
                                gr_meshComm

  use Driver_interface,  ONLY : Driver_abortFlash

  use Grid_interface,    ONLY : Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getBlkBoundBox,    &
                                Grid_getDeltas,         &
                                Grid_getLeafIterator,   &
                                Grid_releaseLeafIterator,&
                                Grid_getCellCoords

  use gr_mpoleData,      ONLY : gr_mpoleDomainXmin,     &
                                gr_mpoleDomainYmin,     &
                                gr_mpoleDomainZmin,     &
                                gr_mpoleDomainXmax,     &
                                gr_mpoleDomainYmax,     &
                                gr_mpoleDomainZmax,     &
                                gr_mpoleDrInnerZone,    &
                                gr_mpoleDrInnerZoneInv, &
                                gr_mpoleXcenter,        &
                                gr_mpoleYcenter,        &
                                gr_mpoleZcenter,        &
                                gr_mpoleTotalMass,      &
                                gr_mpoleXdens2CoM,      &
                                gr_mpoleYdens2CoM,      &
                                gr_mpoleZdens2CoM,      &
                                gr_mpoleXcenterOfMass,  &
                                gr_mpoleYcenterOfMass,  &
                                gr_mpoleZcenterOfMass

  use block_metadata,    ONLY : block_metadata_t
  use leaf_iterator,     ONLY : leaf_iterator_t
  
  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  include "Flash_mpi.h"
  
  integer, intent (in) :: idensvar

  logical :: domainXmax, domainYmax, domainZmax
  logical :: guardCells
  logical :: insideBlock
  logical :: invokeRecv

  
  
  integer :: error
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: k,kmin,kmax
  integer :: maxEdges
  integer :: messageTag
  integer :: nCellsI, nCellsJ, nCellsK
  integer :: nEdgesI, nEdgesJ, nEdgesK

  integer :: locate      (1:3)
  integer :: status      (MPI_STATUS_SIZE)
  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  

  real    :: bndBoxILow
  real    :: bndBoxJLow
  real    :: bndBoxKLow
  real    :: cellDensity, cellMass, cellMassDensity, cellVolume
  real    :: DeltaI, DeltaJ, DeltaK
  real    :: DeltaIHalf, DeltaJHalf, DeltaKHalf
  real    :: localMsum, localMDsum, localMDXsum, localMDYsum, localMDZsum
  real    :: localMXsum, localMYsum, localMZsum
  real    :: totalMassDensityInv
  real    :: x, y, z

  real    :: delta     (1:MDIM)
  real    :: localData (1:8)
  real    :: totalData (1:8)
  real    :: bndBox    (LOW:HIGH,1:MDIM)

  real, allocatable :: shifts   (:,:)
  real, pointer     :: solnData (:,:,:,:)

  integer :: lev
  type(block_metadata_t) :: block
  type(leaf_iterator_t) :: itor
  !
!
!     ...Sum quantities over all locally held leaf blocks.
!
!
  localData(1)   = ZERO
  localData(2)  = ZERO
  localData(3) = ZERO
  localData(4) = ZERO
  localData(5) = ZERO
  localData(6)  = ZERO
  localData(7)  = ZERO
  localData(8)  = ZERO

  
  call Grid_getLeafIterator(itor)
  do while(itor%is_valid())
     call itor%blkMetaData(block)
     lev=block%level
     
     blkLimits=block%limits
     
     call Grid_getBlkBoundBox     (block, bndBox)
     call Grid_getDeltas          (lev, delta)
     call Grid_getBlkPtr          (block, solnData)
     imin = blkLimits (LOW, IAXIS)
     jmin = blkLimits (LOW, JAXIS)
     kmin = blkLimits (LOW, KAXIS)  
     imax = blkLimits (HIGH,IAXIS)
     jmax = blkLimits (HIGH,JAXIS)
     kmax = blkLimits (HIGH,KAXIS)

     DeltaI = delta (IAXIS)
     DeltaJ = delta (JAXIS)
     DeltaK = delta (KAXIS)

     DeltaIHalf = DeltaI * HALF
     DeltaJHalf = DeltaJ * HALF
     DeltaKHalf = DeltaK * HALF

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
     bndBoxKLow = bndBox (LOW,KAXIS)

     cellVolume = DeltaI * DeltaJ * DeltaK
     
     z = bndBoxKLow + DeltaKHalf
     do k = kmin,kmax
        y = bndBoxJLow + DeltaJHalf
        do j = jmin,jmax
           x = bndBoxILow + DeltaIHalf
           do i = imin,imax

              cellDensity     = solnData (idensvar,i,j,k)
              cellMass        = cellDensity * cellVolume
              cellMassDensity = cellMass * cellDensity

              localData(1)   = localData(1)   + cellMass
              localData(2)  = localData(2)  + cellMassDensity
              localData(3) = localData(3) + cellMassDensity * x
              localData(4) = localData(4) + cellMassDensity * y
              localData(5) = localData(5) + cellMassDensity * z
              localData(6)  = localData(6)  + cellMass * x
              localData(7)  = localData(7)  + cellMass * y
              localData(8)  = localData(8)  + cellMass * z
              x = x + DeltaI
           end do
           y = y + DeltaJ
        end do
        z = z + DeltaK
     end do

     call Grid_releaseBlkPtr (block, solnData)
     call itor%next()
  end do
  call Grid_releaseLeafIterator(itor)
!
!
!     ...Prepare for a one-time all reduce call.
!
!

  !
!
!     ...Calculate the total sums and give a copy to each processor.
!
  !
  
  call  MPI_AllReduce (localData,   &
                       totalData,   &
                       8,           &
                       FLASH_REAL,  & 
                       MPI_Sum,     &
                       gr_meshComm, &
                       error        )
!
!
!     ...Analyze total mass obtained. If nonsense, abort.
!
!
  gr_mpoleTotalMass = totalData (1)
     
  if (abs (gr_mpoleTotalMass) < tiny (gr_mpoleTotalMass)) then
      call Driver_abortFlash ('[gr_mpoleCen3Dcartesian] ERROR:  gr_mpoleTotalMass <= 0')
  end if
!
!
!     ...Calculate center of multipole expansion coordinates.
!
!
  totalMassDensityInv = ONE / totalData (2)

  gr_mpoleXcenter = totalData (3) * totalMassDensityInv
  gr_mpoleYcenter = totalData (4) * totalMassDensityInv
  gr_mpoleZcenter = totalData (5) * totalMassDensityInv
  gr_mpoleXdens2CoM = gr_mpoleXcenter
  gr_mpoleYdens2CoM = gr_mpoleYcenter
  gr_mpoleZdens2CoM = gr_mpoleZcenter
  gr_mpoleXcenterOfMass = totalData (6) / totalData (1)
  gr_mpoleYcenterOfMass = totalData (7) / totalData (1)
  gr_mpoleZcenterOfMass = totalData (8) / totalData (1)
!
!
!     ...Find the local blockID to which the center of multipole expansion
!        belongs and place the center on the nearest cell corner. Also at
!        this point we determine the inner zone atomic length, since the
!        inner zone is defined around the center of multipole expansion.
!        Whatever processor is doing the relevant calculation sends its
!        final data (updated center of multipole expansion and inner zone
!        atomic length) to the master, which then broadcasts the info.
!
!
  messageTag = 1
  invokeRecv = .true.

  call Grid_getLeafIterator(itor)
  do while(itor%is_valid())
     call itor%blkMetaData(block)
     
     call Grid_getBlkBoundBox (block, bndBox)
     
     insideBlock =       (gr_mpoleXcenter >= bndBox (LOW ,IAXIS)) &
          .and. (gr_mpoleYcenter >= bndBox (LOW ,JAXIS)) &
          .and. (gr_mpoleZcenter >= bndBox (LOW ,KAXIS)) &
          .and. (gr_mpoleXcenter <  bndBox (HIGH,IAXIS)) &
          .and. (gr_mpoleYcenter <  bndBox (HIGH,JAXIS)) &
          .and. (gr_mpoleZcenter <  bndBox (HIGH,KAXIS))
     
     domainXmax  =       (gr_mpoleXcenter == bndBox (HIGH,IAXIS)) &    ! include (however unlikely) the
          .and. (gr_mpoleXcenter == gr_mpoleDomainXmax)       ! missing X part of the domain
     domainYmax  =       (gr_mpoleYcenter == bndBox (HIGH,JAXIS)) &    ! include (however unlikely) the
          .and. (gr_mpoleYcenter == gr_mpoleDomainYmax)       ! missing Y part of the domain
     domainZmax  =       (gr_mpoleZcenter == bndBox (HIGH,KAXIS)) &    ! include (however unlikely) the
          .and. (gr_mpoleZcenter == gr_mpoleDomainZmax)       ! missing Z part of the domain
     
     insideBlock = insideBlock .or. domainXmax .or. domainYmax .or. domainZmax
     
     if (insideBlock) then
        lev=block%level
        call Grid_getDeltas          (lev, delta)
        blkLimits=block%limits
        
        DeltaI = delta (IAXIS)
        DeltaJ = delta (JAXIS)
        DeltaK = delta (KAXIS)
        
        gr_mpoleDrInnerZone = HALF * (DeltaI * DeltaJ * DeltaK) ** (ONE / THREE)
        imin = blkLimits (LOW, IAXIS)
        jmin = blkLimits (LOW, JAXIS)
        kmin = blkLimits (LOW, KAXIS)  
        imax = blkLimits (HIGH,IAXIS)
        jmax = blkLimits (HIGH,JAXIS)
        kmax = blkLimits (HIGH,KAXIS)
        
        nCellsI = imax - imin + 1
        nCellsJ = jmax - jmin + 1
        nCellsK = kmax - kmin + 1
        
        nEdgesI = nCellsI + 1
        nEdgesJ = nCellsJ + 1
        nEdgesK = nCellsK + 1
        
        maxEdges = max (nEdgesI, nEdgesJ, nEdgesK)
        
        allocate (shifts (1:maxEdges,3))
        
        guardCells = .false.
        
        call Grid_getCellCoords (IAXIS, block, FACES, guardCells, shifts (1:nEdgesI,1), nEdgesI)
        call Grid_getCellCoords (JAXIS, block, FACES, guardCells, shifts (1:nEdgesJ,2), nEdgesJ)
        call Grid_getCellCoords (KAXIS, block, FACES, guardCells, shifts (1:nEdgesK,3), nEdgesK)
        
        shifts (1:nEdgesI,1) = shifts (1:nEdgesI,1) - gr_mpoleXcenter
        shifts (1:nEdgesJ,2) = shifts (1:nEdgesJ,2) - gr_mpoleYcenter
        shifts (1:nEdgesK,3) = shifts (1:nEdgesK,3) - gr_mpoleZcenter
        
        locate (1) = minloc (abs (shifts (1:nEdgesI,1)), dim = 1)
        locate (2) = minloc (abs (shifts (1:nEdgesJ,2)), dim = 1)
        locate (3) = minloc (abs (shifts (1:nEdgesK,3)), dim = 1)
        
        gr_mpoleXcenter = gr_mpoleXcenter + shifts (locate (1),1)  ! move to nearest x edge
        gr_mpoleYcenter = gr_mpoleYcenter + shifts (locate (2),2)  ! move to nearest y edge
        gr_mpoleZcenter = gr_mpoleZcenter + shifts (locate (3),3)  ! move to nearest z edge
        deallocate (shifts)
        
        localData (1) = gr_mpoleDrInnerZone
        localData (2) = gr_mpoleXcenter
        localData (3) = gr_mpoleYcenter
        localData (4) = gr_mpoleZcenter
        
        if (gr_meshMe /= MASTER_PE) then
           
           call MPI_Send (localData,    &
                4,            &
                FLASH_REAL,   &
                MASTER_PE,    &
                messageTag,   &
                gr_meshComm,  &
                error         )
        else
           invokeRecv = .false.
        end if
        
        exit
        
     end if
     call itor%next()
  end do
  call Grid_releaseLeafIterator(itor)
  if ((gr_meshMe == MASTER_PE) .and. invokeRecv) then
     
     call MPI_Recv (localData,      &
          4,              &
          FLASH_REAL,     &
          MPI_ANY_SOURCE, &
          messageTag,     &
          gr_meshComm,    &
          status,         &
          error           )
  end if
  !
  !
  !     ...At this point, the master has all the info. Broadcast and update all
  !        other processors.
  !
  !
  call MPI_Bcast (localData,   &
       4,           &
       FLASH_REAL,  &
       MASTER_PE,   &
       gr_meshComm, &
       error        )
  
  gr_mpoleDrInnerZone    = localData (1)
  gr_mpoleDrInnerZoneInv = ONE / gr_mpoleDrInnerZone
  gr_mpoleXcenter        = localData (2)
  gr_mpoleYcenter        = localData (3)
  gr_mpoleZcenter        = localData (4)
  !
  !
  !     ...Ready!
  !
  !
  return
end subroutine gr_mpoleCen3Dcartesian
