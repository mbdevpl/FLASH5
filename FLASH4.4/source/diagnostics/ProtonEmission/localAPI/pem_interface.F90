!!****ih* source/diagnostics/ProtonEmission/localAPI/pem_interface
!!
!! NAME
!!
!!  pem_interface
!!
!! SYNOPSIS
!!
!!   use pem_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Proton Emission unit that defines its
!!  private interfaces.
!!
!!***

Module pem_interface

  interface
     logical function pem_activeEmissionRegion (x,y,z)
       real, intent (in) :: x,y,z
     end function pem_activeEmissionRegion
  end interface

  interface
     real function pem_avReactivityDdpT (TkeV)
       real, intent (in) :: TkeV
     end function pem_avReactivityDdpT
  end interface

  interface
     real function pem_avReactivityHe3dpHe4 (TkeV)
       real, intent (in) :: TkeV
     end function pem_avReactivityHe3dpHe4
  end interface

  interface
     subroutine pem_blockData3DRec (iminBlock, imaxBlock,            &
                                    jminBlock, jmaxBlock,            &
                                    kminBlock, kmaxBlock,            &
                                    deltaInvX, deltaInvY, deltaInvZ, &
                                    blockData                        )

       integer, intent (in) :: iminBlock, imaxBlock
       integer, intent (in) :: jminBlock, jmaxBlock
       integer, intent (in) :: kminBlock, kmaxBlock
       real,    intent (in) :: deltaInvX, deltaInvY, deltaInvZ
       real,    intent (in) :: blockData (:,:,:,:)
     end subroutine pem_blockData3DRec
  end interface

  interface
     subroutine pem_closeDetectorFiles () 
     end subroutine pem_closeDetectorFiles
  end interface

  interface
     subroutine pem_closeEmissionProfileFile () 
     end subroutine pem_closeEmissionProfileFile
  end interface

  interface
     subroutine pem_createBlockProtons3DRec (blockID,              &
                                             iminBlock, imaxBlock, &
                                             jminBlock, jmaxBlock, &
                                             kminBlock, kmaxBlock, &
                                             timeStep,             &
                                             cellVolume,           &
                                             countOnly,            &
                                             blockData,            &
                                             blockIncomplete       )
       integer, intent (in)    :: blockID
       integer, intent (in)    :: iminBlock, imaxBlock
       integer, intent (in)    :: jminBlock, jmaxBlock
       integer, intent (in)    :: kminBlock, kmaxBlock
       real,    intent (in)    :: timeStep
       real,    intent (in)    :: cellVolume
       logical, intent (in)    :: countOnly
       real,    intent (in)    :: blockData (:,:,:,:)
       logical, intent (inout) :: blockIncomplete
     end subroutine pem_createBlockProtons3DRec
  end interface

  interface
     subroutine pem_createProtons3DRec (blockCount,                  &
                                        blockList,                   &
                                        timeStep,                    &
                                        countOnly,                   &
                                                   emissionInitiate, &
                                                   emissionIncomplete)

       integer, intent (in)    :: blockCount
       integer, intent (in)    :: blockList (1:blockCount)
       real,    intent (in)    :: timeStep
       logical, intent (in)    :: countOnly
       logical, intent (inout) :: emissionInitiate
       logical, intent (inout) :: emissionIncomplete
     end subroutine pem_createProtons3DRec
  end interface

  interface
     subroutine pem_createProtons (blockCount,                  &
                                   blockList,                   &
                                   timeStep,                    &
                                   countOnly,                   &
                                              emissionInitiate, &
                                              emissionIncomplete)

       integer, intent (in)    :: blockCount
       integer, intent (in)    :: blockList (1:blockCount)
       real,    intent (in)    :: timeStep
       logical, intent (in)    :: countOnly
       logical, intent (inout) :: emissionInitiate
       logical, intent (inout) :: emissionIncomplete
     end subroutine pem_createProtons
  end interface

  interface
     subroutine pem_createProtonTags () 
     end subroutine pem_createProtonTags
  end interface

  interface
     subroutine pem_detectorsCheck3DRec () 
     end subroutine pem_detectorsCheck3DRec
  end interface

  interface
     subroutine pem_detectorsCheck () 
     end subroutine pem_detectorsCheck
  end interface

  interface
     subroutine pem_detectorsInfo3DRec () 
     end subroutine pem_detectorsInfo3DRec
  end interface

  interface
     subroutine pem_detectorsInfo () 
     end subroutine pem_detectorsInfo
  end interface

  interface
     subroutine pem_flushScreenProtons2Disk () 
     end subroutine pem_flushScreenProtons2Disk
  end interface

  interface
     function pem_maxConfinement3DRec (nc,y)
       integer, intent (in) :: nc
       real,    intent (in) :: y (:)
       real                 :: pem_maxConfinement3DRec (1:nc)
     end function pem_maxConfinement3DRec
  end interface

  interface
     function pem_minConfinement3DRec (nc,y)
       integer, intent (in) :: nc
       real,    intent (in) :: y (:)
       real                 :: pem_minConfinement3DRec (1:nc)
     end function pem_minConfinement3DRec
  end interface

  interface
     subroutine pem_openDetectorFiles () 
     end subroutine pem_openDetectorFiles
  end interface

  interface
     subroutine pem_openEmissionProfileFile () 
     end subroutine pem_openEmissionProfileFile
  end interface

  interface
     subroutine pem_printBlockVariable (blockID, variable, fileUnit)
       integer, intent (in) :: blockID
       integer, intent (in) :: variable
       integer, intent (in) :: fileUnit
     end subroutine pem_printBlockVariable
  end interface

  interface
     subroutine pem_printDetectorsData () 
     end subroutine pem_printDetectorsData
  end interface

  interface
     subroutine pem_printEmissionStamp (timeStep, timeSimulation)
       real, intent (in) :: timeStep
       real, intent (in) :: timeSimulation
     end subroutine pem_printEmissionStamp
  end interface

  interface
     subroutine pem_printMainData () 
     end subroutine pem_printMainData
  end interface

  interface
     subroutine pem_printMatrix (fileUnit,                   &
                                 title,                      &
                                 rowMinMatrix, rowMaxMatrix, &
                                 colMinMatrix, colMaxMatrix, &
                                 rowMinPrint , rowMaxPrint,  &
                                 colMinPrint , colMaxPrint,  &
                                 matrix                      )

      integer,           intent (in) :: fileUnit
      character (len=*), intent (in) :: title
      integer,           intent (in) :: rowMinMatrix, rowMaxMatrix
      integer,           intent (in) :: colMinMatrix, colMaxMatrix
      integer,           intent (in) :: rowMinPrint , rowMaxPrint
      integer,           intent (in) :: colMinPrint , colMaxPrint
      real,              intent (in) :: matrix (rowMinMatrix : rowMaxMatrix , colMinMatrix : colMaxMatrix)
     end subroutine pem_printMatrix
  end interface

  interface
     subroutine pem_printProtonsData (fileLabel, processorID) 
       character (len=*), intent (in) :: fileLabel
       integer,           intent (in) :: processorID
     end subroutine pem_printProtonsData
  end interface

  interface
     subroutine pem_printProtonSourcesData () 
     end subroutine pem_printProtonSourcesData
  end interface

  interface
     subroutine pem_protonsBlockIDInfo () 
     end subroutine pem_protonsBlockIDInfo
  end interface

  interface
     real function pem_protonSpeed (EMeV)
       real, intent (in) :: EMeV
     end function pem_protonSpeed
  end interface

  interface
     subroutine pem_recordProtonOnScreen (px, py, pz, vx, vy, vz)
       real,    intent (in) :: px, py, pz
       real,    intent (in) :: vx, vy, vz
     end subroutine pem_recordProtonOnScreen
  end interface

  interface
     subroutine pem_rotationMatrix3DaxisZ (Vx,Vy,Vz,  R11,R12,R13,R21,R22,R23,R31,R32,R33)
       real, intent (in)  :: Vx, Vy, Vz
       real, intent (out) :: R11,R12,R13
       real, intent (out) :: R21,R22,R23
       real, intent (out) :: R31,R32,R33
     end subroutine pem_rotationMatrix3DaxisZ
  end interface

  interface
     subroutine pem_setupDetectors () 
     end subroutine pem_setupDetectors
  end interface

  interface
     subroutine pem_setupEmissionBoxes () 
     end subroutine pem_setupEmissionBoxes
  end interface

  interface
     subroutine pem_setupEmissionBuckets () 
     end subroutine pem_setupEmissionBuckets
  end interface

  interface
     subroutine pem_setupProtons () 
     end subroutine pem_setupProtons
  end interface

  interface
     subroutine pem_setupProtonSources () 
     end subroutine pem_setupProtonSources
  end interface

  interface
     subroutine pem_setupScreenProtons () 
     end subroutine pem_setupScreenProtons
  end interface

  interface
     subroutine pem_statisticalFinalize ()
     end subroutine pem_statisticalFinalize
  end interface

  interface
     subroutine pem_statisticalInit ()
     end subroutine pem_statisticalInit
  end interface

  interface
     subroutine pem_statisticalSetSeed (seed)
       integer, intent (in) :: seed
     end subroutine pem_statisticalSetSeed
  end interface

  interface
     subroutine pem_statisticalXYZunitSphconeZ (cosAlpha, sinAlpha,       &
                                                arraySize,                &
                                                nWanted,                  &
                                                               nReturned, &
                                                               xSphcone,  &
                                                               ySphcone,  &
                                                               zSphcone   )

       real,    intent (in)  :: cosAlpha, sinAlpha
       integer, intent (in)  :: arraySize
       integer, intent (in)  :: nWanted
       integer, intent (out) :: nReturned
       real,    intent (out) :: xSphcone (1:arraySize)
       real,    intent (out) :: ySphcone (1:arraySize)
       real,    intent (out) :: zSphcone (1:arraySize)
     end subroutine pem_statisticalXYZunitSphconeZ
  end interface

  interface
     subroutine pem_statisticalXYZunitSphere (arraySize,          &
                                                         nSphere, &
                                                         xSphere, &
                                                         ySphere, &
                                                         zSphere  )

       integer, intent (in)  :: arraySize
       integer, intent (out) :: nSphere
       real,    intent (out) :: xSphere (1:arraySize)
       real,    intent (out) :: ySphere (1:arraySize)
       real,    intent (out) :: zSphere (1:arraySize)
     end subroutine pem_statisticalXYZunitSphere
  end interface

  interface
     real function pem_time2FacesParabolicPath1D (pos, vel, acc, minFace, maxFace, noFaceTime)
       real, intent (in) :: pos, vel, acc
       real, intent (in) :: minFace, maxFace
       real, intent (in) :: noFaceTime
     end function pem_time2FacesParabolicPath1D
  end interface

  interface
     subroutine pem_traceBlockProtons3DRec (protonFirst,  protonLast,          &
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
     end subroutine pem_traceBlockProtons3DRec
  end interface

  interface
     function pem_traceODEfunction3DRec (t,y)
       real, intent (in) :: t
       real, intent (in) :: y (:)
       real :: pem_traceODEfunction3DRec (1:size (y))
     end function pem_traceODEfunction3DRec
  end interface

  interface
     subroutine pem_traceProtons3DRec () 
     end subroutine pem_traceProtons3DRec
  end interface

  interface
     subroutine pem_traceProtons () 
     end subroutine pem_traceProtons
  end interface

  interface
     subroutine pem_updateProtons (doMove) 
       logical, intent (in) :: doMove
     end subroutine pem_updateProtons
  end interface

  interface
     subroutine pem_vectorOrthoBasis3D (Vx,Vy,Vz,  ix,iy,iz,jx,jy,jz,kx,ky,kz)
       real, intent (in)  :: Vx, Vy, Vz
       real, intent (out) :: ix, iy, iz
       real, intent (out) :: jx, jy, jz
       real, intent (out) :: kx, ky, kz
     end subroutine pem_vectorOrthoBasis3D
  end interface

  interface
     subroutine pem_write2DetectorFile (detector, recordCount, bucketCount, bucket)
       integer, intent (in) :: detector
       integer, intent (in) :: recordCount
       integer, intent (in) :: bucketCount
       real,    intent (in) :: bucket (1:recordCount,1:bucketCount)
     end subroutine pem_write2DetectorFile
  end interface

end Module pem_interface
