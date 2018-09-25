!!****if* source/physics/Eos/unitTest/eos_initTest
!!
!! NAME
!!
!!  eos_initTest
!!
!!
!! SYNOPSIS
!!
!!  call eos_initTest()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize Eos unitTest scope variables which are typically from runtime parameters.
!!  The part of the unitTest located in the Simulation unit tree should arrange for this
!!  initialization routine to be called.
!!
!! ARGUMENTS
!!
!!
!!  none
!!  
!!
!! PARAMETERS
!!
!!
!!    eos_testPresMode [STRING]
!!    eos_testEintMode [STRING]
!!    eos_testTempMode [STRING]
!!***

subroutine eos_initTest()

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
                                         RuntimeParameters_mapStrToInt
  use eos_testData, ONLY: eos_testPresModeStr, &
                          eos_testEintModeStr, &
                          eos_testTempModeStr, &
                          eos_testPresMode, &
                          eos_testEintMode, &
                          eos_testTempMode
  use eos_testData, ONLY: eos_testTolerance
  use eos_testData, ONLY: eos_test1allB, &
                          eos_test2allB, &
                          eos_test3allB, &
                          eos_test4allB

  implicit none

  call RuntimeParameters_get ("eos_testPresMode", eos_testPresModeStr)
  call RuntimeParameters_mapStrToInt(eos_testPresModeStr, eos_testPresMode)
  call RuntimeParameters_get ("eos_testEintMode", eos_testEintModeStr)
  call RuntimeParameters_mapStrToInt(eos_testEintModeStr, eos_testEintMode)
  call RuntimeParameters_get ("eos_testTempMode", eos_testTempModeStr)
  call RuntimeParameters_mapStrToInt(eos_testTempModeStr, eos_testTempMode)

  call RuntimeParameters_get ("eos_testTolerance", eos_testTolerance)

  eos_test1allB = .TRUE.
  eos_test2allB = .TRUE.
  eos_test3allB = .TRUE.
  eos_test4allB = .TRUE.

end subroutine eos_initTest
