!!****ih* source/diagnostics/ThomsonScattering/localAPI/thsc_interface
!!
!! NAME
!!
!!  thsc_interface
!!
!! SYNOPSIS
!!
!!   use thsc_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Thomson Scattering unit that defines its
!!  private interfaces.
!!
!!***

#include "constants.h"

Module thsc_interface

  use thsc_constInterface ! for integration results: (optical depth, Faraday rotation angle, etc.)
  use thsc_setupDetectedRays_mod


  interface
     subroutine thsc_beamsCheck2DCyl3D () 
     end subroutine thsc_beamsCheck2DCyl3D
  end interface

  interface
     subroutine thsc_beamsCheck3DRec () 
     end subroutine thsc_beamsCheck3DRec
  end interface

  interface
     subroutine thsc_beamsCheck () 
     end subroutine thsc_beamsCheck
  end interface

  interface
     subroutine thsc_beamsInfo3DRec () 
     end subroutine thsc_beamsInfo3DRec
  end interface

  interface
     subroutine thsc_beamsInfo () 
     end subroutine thsc_beamsInfo
  end interface

  interface
     subroutine thsc_blockData3DRec (iminBlock, imaxBlock,            &
                                   jminBlock, jmaxBlock,            &
                                   kminBlock, kmaxBlock,            &
                                   deltaInvX, deltaInvY, deltaInvZ, &
                                   blockData                        )

       integer, intent (in) :: iminBlock, imaxBlock
       integer, intent (in) :: jminBlock, jmaxBlock
       integer, intent (in) :: kminBlock, kmaxBlock
       real,    intent (in) :: deltaInvX, deltaInvY, deltaInvZ
       real,    intent (in) :: blockData (:,:,:,:)
     end subroutine thsc_blockData3DRec
  end interface

  interface
     subroutine thsc_calculateCellData (coordMin, coordMax, indexMin, indexMax, delta, coords)
       real,    intent (in)  :: coordMin, coordMax
       integer, intent (in)  :: indexMin, indexMax
       real,    intent (out) :: delta
       real,    intent (out) :: coords (indexMin:indexMax+1)
     end subroutine thsc_calculateCellData
  end interface

  interface
     subroutine thsc_capsuleGrainIndices2xyz (i,j,k,S,x,y,z)
       integer, intent (in)  :: i,j,k
       real,    intent (in)  :: S
       real,    intent (out) :: x,y,z
     end subroutine thsc_capsuleGrainIndices2xyz
  end interface

  interface
     subroutine thsc_capsuleNextGrainIndices (L,i,j,k,valid)
       integer, intent (in)    :: L
       integer, intent (inout) :: i,j,k
       logical, intent (out)   :: valid
     end subroutine thsc_capsuleNextGrainIndices
  end interface

  interface
     integer function thsc_capsuleTotalGrainCount (L)
       integer, intent (in) :: L
     end function thsc_capsuleTotalGrainCount
  end interface

  interface thsc_checkCriticalFreq
     subroutine thsc_checkCriticalFreq(nele, omgL, forbidden)
       implicit none
       real   , intent(in)  :: nele
       real   , intent(in)  :: omgL
       logical, intent(OUT) :: forbidden
     end subroutine thsc_checkCriticalFreq

     subroutine thsc_checkCriticalFreqFullState(solnVec, omgL, forbidden)
       implicit none
       real   , intent(in)  :: solnVec(:)
       real   , intent(in)  :: omgL
       logical, intent(OUT) :: forbidden
     end subroutine thsc_checkCriticalFreqFullState
  end interface

  interface
     subroutine thsc_closeDetectorFiles () 
     end subroutine thsc_closeDetectorFiles
  end interface

  interface
     subroutine thsc_computeBeamPower (timeStep, timeSimulation, pulseNumber, beamPower)
       implicit none
       real,    intent (in)  :: timeStep
       real,    intent (in)  :: timeSimulation
       integer, intent (in)  :: pulseNumber
       real,    intent (out) :: beamPower
     end subroutine thsc_computeBeamPower
  end interface

  interface
     subroutine thsc_createProtons3DRec (blockCount, blockList, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeSimulation
     end subroutine thsc_createProtons3DRec
  end interface

  interface
     subroutine thsc_createProtons (blockCount, blockList, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeSimulation
     end subroutine thsc_createProtons
  end interface

  interface
     subroutine thsc_detectorsCheck3DRec () 
     end subroutine thsc_detectorsCheck3DRec
  end interface

  interface
     subroutine thsc_detectorsCheck () 
     end subroutine thsc_detectorsCheck
  end interface

  interface
     subroutine thsc_detectorsInfo3DRec () 
     end subroutine thsc_detectorsInfo3DRec
  end interface

  interface
     subroutine thsc_detectorsInfo () 
     end subroutine thsc_detectorsInfo
  end interface

  interface
     subroutine thsc_flushDetectorData2Disk () 
     end subroutine thsc_flushDetectorData2Disk
  end interface

  interface
     subroutine thsc_initializeBeams () 
     end subroutine thsc_initializeBeams
  end interface

  interface
     subroutine thsc_initializeDetectors () 
     end subroutine thsc_initializeDetectors
  end interface

  interface
     function thsc_maxConfinement3DRec (nc,y)
       integer, intent (in) :: nc
       real,    intent (in) :: y (:)
       real                 :: thsc_maxConfinement3DRec (1:nc)
     end function thsc_maxConfinement3DRec
  end interface

  interface
     function thsc_minConfinement3DRec (nc,y)
       integer, intent (in) :: nc
       real,    intent (in) :: y (:)
       real                 :: thsc_minConfinement3DRec (1:nc)
     end function thsc_minConfinement3DRec
  end interface

  interface
     subroutine thsc_openDetectorFiles (timeSimulation) 
       real, intent (in) :: timeSimulation
     end subroutine thsc_openDetectorFiles
  end interface

  interface
     real function thsc_parabolicPathLength3D (vx,vy,vz,ax,ay,az,T)
       real, intent (in) :: vx,vy,vz
       real, intent (in) :: ax,ay,az
       real, intent (in) :: T
     end function thsc_parabolicPathLength3D
  end interface

  interface
     subroutine thsc_printBeamsData () 
     end subroutine thsc_printBeamsData
  end interface

  interface
     subroutine thsc_printBlockVariable (blockID, variable, fileUnit)
       integer, intent (in) :: blockID
       integer, intent (in) :: variable
       integer, intent (in) :: fileUnit
     end subroutine thsc_printBlockVariable
  end interface

  interface
     subroutine thsc_printDetectorsData () 
     end subroutine thsc_printDetectorsData
  end interface

  interface
     subroutine thsc_printMainData () 
     end subroutine thsc_printMainData
  end interface

  interface
     subroutine thsc_printPulsesData () 
     end subroutine thsc_printPulsesData
  end interface

  interface
     subroutine thsc_printMatrix (fileUnit,                   &
                                title,                      &
                                rowMinMatrix, rowMaxMatrix, &
                                colMinMatrix, colMaxMatrix, &
                                rowMinPrint , rowMaxPrint,  &
                                colMinPrint , colMaxPrint,  &
                                matrix                      )
       implicit none
       integer,           intent (in) :: fileUnit
       character (len=*), intent (in) :: title
       integer,           intent (in) :: rowMinMatrix, rowMaxMatrix
       integer,           intent (in) :: colMinMatrix, colMaxMatrix
       integer,           intent (in) :: rowMinPrint , rowMaxPrint
       integer,           intent (in) :: colMinPrint , colMaxPrint
       real,              intent (in) :: matrix (rowMinMatrix : rowMaxMatrix , colMinMatrix : colMaxMatrix)
     end subroutine thsc_printMatrix
  end interface

  interface
     subroutine thsc_printProtonsData (fileLabel, processorID) 
       character (len=*), intent (in) :: fileLabel
       integer,           intent (in) :: processorID
     end subroutine thsc_printProtonsData
  end interface

  interface
     subroutine thsc_protonsBlockIDInfo () 
     end subroutine thsc_protonsBlockIDInfo
  end interface

  INTERFACE
     SUBROUTINE THSC_PRUNE1(BLOCKCOUNT,BLOCKLIST)
       INTEGER, INTENT(IN) :: BLOCKCOUNT
       INTEGER, INTENT(IN) :: BLOCKLIST(1:BLOCKCOUNT)
     END SUBROUTINE THSC_PRUNE1
     SUBROUTINE THSC_PRUNE2(BLOCKCOUNT,BLOCKLIST)
       INTEGER, INTENT(IN) :: BLOCKCOUNT
       INTEGER, INTENT(IN) :: BLOCKLIST(1:BLOCKCOUNT)
     END SUBROUTINE THSC_PRUNE2
  END INTERFACE

  interface
     subroutine thsc_recordRayOnDetector(solnVec, &
          power, &
          vpar, vperp, &
          sa, &
          ud, gamma, &
!!$          px, py, pz, vx, vy, vz, Jv, Kx, Ky, Kz, &
          scatterCoords, dv, &
          beam, detector,    &
          xl, yl, zl)
!!$       real,    intent (in) :: px, py, pz
!!$       real,    intent (in) :: vx, vy, vz
!!$       real,    intent (in) :: Jv, Kx, Ky, Kz
       real,    intent(in)  :: solnVec(:)
       real,    intent (in) :: power
       real,    intent (in) :: vpar, vperp, ud
       real,    intent (in) :: sa, gamma
       real,    intent (in) :: scatterCoords(MDIM)
       real,    intent (in) :: dv
       integer, intent (in) :: detector, beam
       real,    intent (in) :: xl,yl,zl
     end subroutine thsc_recordRayOnDetector
  end interface

  interface
     subroutine thsc_setupBeams () 
     end subroutine thsc_setupBeams
  end interface

  interface
     subroutine thsc_setupPulses () 
     end subroutine thsc_setupPulses
  end interface

  interface
     subroutine thsc_setupDetectors () 
     end subroutine thsc_setupDetectors
  end interface

  interface
     subroutine thsc_setupRays () 
     end subroutine thsc_setupRays
  end interface

!!$  interface
!!$     subroutine thsc_setupDetectedRays () 
!!$     end subroutine thsc_setupDetectedRays
!!$  end interface

  interface
     real function thsc_time2FacesParabolicPath1D (pos, vel, acc, minFace, maxFace, noFaceTime)
       real, intent (in) :: pos, vel, acc
       real, intent (in) :: minFace, maxFace
       real, intent (in) :: noFaceTime
     end function thsc_time2FacesParabolicPath1D
  end interface

  interface
     subroutine thsc_rayLoop (blockCount, blockList, simTime)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: simTime
     end subroutine thsc_rayLoop
  end interface

  interface
     subroutine thsc_rayGetPowerIn(xs,ys,zs,xc,yc,zc,kiHatX,kiHatY,kiHatZ,ibeam,powerIn,inBeam, &
                                   inMeasureRegion,&
                                   xl,yl,zl)
       implicit none
       real, intent(in)  :: xs, ys, zs
       real, intent(INOUT) :: xc, yc, zc
       real, intent(OUT) :: kiHatX,kiHatY,kiHatZ
       integer, intent(in) :: ibeam
       real, intent(OUT) :: powerIn
       logical, intent(OUT) :: inBeam
       logical, intent(OUT) :: inMeasureRegion
       real, intent(OUT),OPTIONAL :: xl,yl,zl
     end subroutine thsc_rayGetPowerIn
  end interface

  interface
     subroutine thsc_integrate (blockCount, blockList, leg)
       implicit none
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       integer, intent (in) :: leg
     end subroutine thsc_integrate
  end interface

  interface
     subroutine thsc_traceToDetector(solnVec,xs,ys,zs,xc,yc,zc,kiHatX,kiHatY,kiHatZ, dv, power, &
          ibeam,detector,&
          inDetectorCone, xl,yl,zl)
       implicit none
       real, intent(INOUT)  :: solnVec(:)
       real, intent(in)  :: xs, ys, zs
       real, intent(in)  :: xc, yc, zc
       real, intent(in)  :: kiHatX,kiHatY,kiHatZ
       real, intent(in)  :: dv
       real, intent(in)  :: power
       integer, intent(in) :: ibeam, detector
       logical, intent(OUT) :: inDetectorCone
       real, intent(in)  :: xl,yl,zl
     end subroutine thsc_traceToDetector
  end interface

  interface
     subroutine thsc_traceBlockProtons3DRec (protonFirst,  protonLast,          &
                                           iminBlock, imaxBlock,              &
                                           jminBlock, jmaxBlock,              &
                                           kminBlock, kmaxBlock,              &
                                           xminBlock, xmaxBlock,              &
                                           yminBlock, ymaxBlock,              &
                                           zminBlock, zmaxBlock,              &
                                           deltaX, deltaY, deltaZ,            &
                                           deltaInvX, deltaInvY, deltaInvZ,   &
                                           blockReflectMinX,                  &
                                           blockReflectMaxX,                  &
                                           blockReflectMinY,                  &
                                           blockReflectMaxY,                  &
                                           blockReflectMinZ,                  &
                                           blockReflectMaxZ                   ) 

       integer, intent (in)    :: protonFirst, protonLast   
       integer, intent (in)    :: iminBlock, imaxBlock
       integer, intent (in)    :: jminBlock, jmaxBlock
       integer, intent (in)    :: kminBlock, kmaxBlock
       real,    intent (in)    :: xminBlock, xmaxBlock
       real,    intent (in)    :: yminBlock, ymaxBlock
       real,    intent (in)    :: zminBlock, zmaxBlock
       real,    intent (in)    :: deltaX, deltaY, deltaZ
       real,    intent (in)    :: deltaInvX, deltaInvY, deltaInvZ
       logical, intent (in)    :: blockReflectMinX
       logical, intent (in)    :: blockReflectMaxX
       logical, intent (in)    :: blockReflectMinY
       logical, intent (in)    :: blockReflectMaxY
       logical, intent (in)    :: blockReflectMinZ
       logical, intent (in)    :: blockReflectMaxZ
     end subroutine thsc_traceBlockProtons3DRec
  end interface

  interface
     function thsc_traceODEfunction3DRec (t,y)
       real, intent (in) :: t
       real, intent (in) :: y (:)
       real :: thsc_traceODEfunction3DRec (1:size (y))
     end function thsc_traceODEfunction3DRec
  end interface

  interface
     subroutine thsc_write2DetectorFile (detector, recordCount, bucketCount, bucket)
       integer, intent (in) :: detector
       integer, intent (in) :: recordCount
       integer, intent (in) :: bucketCount
       real,    intent (INOUT) :: bucket (1:recordCount,1:bucketCount)
     end subroutine thsc_write2DetectorFile
  end interface

  interface
     subroutine thsc_write2SpectrumFile (detector, npts, omgs, spectrum)
       integer, intent (in) :: detector
       integer, intent (in) :: npts
       real,    intent (in) :: omgs(1:npts)
       real,    intent (in) :: spectrum(1:npts)
     end subroutine thsc_write2SpectrumFile
  end interface

  interface
     subroutine thsc_ThomsonFormFactor (nFreqs, nAtoms, &
                                      Tele, Tion,     &
                                      Z, A, fract,    &
                                      Nele,           &
                                      Vpar, Vperp,    &
                                      ud,             &
                                      omgs, omgL,     &
                                      sa, dphi, gamma,&
                                      formFactor)
       implicit none
       integer, intent (in) :: nFreqs
       integer, intent (in) :: nAtoms
       real,    intent (in) :: Tele, Tion
       real,    intent (in) :: Z     (1:nAtoms)
       real,    intent (in) :: A     (1:nAtoms)
       real,    intent (in) :: fract (1:nAtoms)
       real,    intent (in) :: Nele
       real,    intent (in) :: Vpar, Vperp
       real,    intent (in) :: ud
       real,    intent (in) :: omgs  (1:nFreqs), omgL
       real,    intent (OUT):: formFactor(1:nFreqs)
       real,    intent (in) :: sa, dphi, gamma
     end subroutine thsc_ThomsonFormFactor
  end interface

  interface
     subroutine thsc_loadExternalData ()
     end subroutine thsc_loadExternalData
  end interface
end Module thsc_interface
