!!****if* source/physics/Eos/unitTest/Eos_unitTest
!! NAME
!!
!!  Eos_unitTest
!! 
!! SYNOPSIS
!!
!!  call Eos_unitTest(integer(IN) :: fileUnit,
!!                    logical(OUT) :: perfect
!!
!! DESCRIPTION
!!
!! This function is the unit test for the Eos unit. It is invoked in
!! the setup unitTest/Eos. The Config file for Eos unit test setup
!! requests a few extra variables in the main grid data structure for
!! Grid scope temporary storage. The Simulation_initBlock of the Eos
!! unit test initializes density in the right place for the DENS_VAR
!! variable (see Flash.h for DENS_VAR, TEMP_VAR etc definitions), and
!! temperature and pressure in the extra storage space CTMP_VAR
!! and CPRS_VAR. The physical quantities at this point are not in
!! thermal equilibrium. 
!!
!! The Eos_unit test starts by copying the initialized
!! temperature into the TEMP_VAR location and calling the
!! Eos_wrapped function with eosMode = MODE_DENS_TEMP, where
!! density and temperature are given and pressure and energy are
!! calculated. Now PRES_VAR and EINT_VAR contain values of pressure
!! and internal energy that are in thermal equilibrium, and the pressure values
!! are not necessarily what was stored in the extra storage space
!! during intialization. 
!! 
!! At this point in time three quantities; temperature,
!! pressure and energy are saved in the extra storage requested by
!! the unitTest/Eos setup, say OTMP_VAR, OPRS_VAR and OENT_VAR. Now
!! the Eos_unitTest function calls Eos_wrapped with eosMode =
!! MODE_DENS_PRES, followed by eosMode= MODE_DENS_EI.  If the
!! newly calculated values of temperature, pressure and energy are
!! the same as those saved in OTMP_VAR, OPRS_VAR and OENT_VAR, then
!! we can conclude that the Eos is working in MODE_DENS_PRES and
!! MODE_DENS_EI modes. However, we still can't say anything about the
!! MODE_DENS_TEMP mode. So we repeat the process by copying CPRS_VAR
!! into PRES_VAR and calling Eos_wrapped with MODE_DENS_PRES. We
!! again save the calculated values in the extra storage and make two
!! more Eos_wrapped calls with the remaining two modes. This time if
!! the new and old values of variables compare, we can conclude that
!! MODE_DENS_TEMP works too, and hence the unit test is successful.
!!
!! A final test calculates the optional derivates by setting the mask
!! argument to Eos equal to true.  The values of these arguments are not
!! tested in any way.  The unitTest simply makes sure that they can be calculated
!! without NaNs or the like.
!!
!!  ARGUMENTS 
!!   
!! 
!!   fileUnit : unit number for file opened by the unitTest/Eos setup
!!              in which to write results of the test
!!
!!   perfect : indicates test ran without error is true.
!!
!!   blockDesc : describes a block on which we are being called.
!!               For informational messages ony, can be omitted.
!!
!!  PARAMETERS
!!
!!  eintSwitch  a rarely used switch which ensures that internal energy calculations 
!!        maintain sufficient precision. Important only if energyTotal is dominated 
!!        by energyKinetic.
!!
!!***

!!REORDER(4): solnData

subroutine Eos_unitTest4(fileUnit, perfect, solnData, blkLimits, blockDesc)

  use Eos_interface, ONLY : Eos_wrapped, Eos

  use block_metadata, ONLY : block_metadata_t
  use Eos_data, ONLY : eos_meshMe, eos_meshNumProcs
  use eos_testData, ONLY: eos_testPresModeStr, &
                          eos_testEintModeStr, &
                          eos_testTempModeStr, &
                          eos_testPresMode, &
                          eos_testEintMode, &
                          eos_testTempMode
  use eos_testData, ONLY: tolerance => eos_testTolerance
  ! logical flags meaning partiall success for all blocks:
  use eos_testData, ONLY: test1allB => eos_test1allB, &
                          test2allB => eos_test2allB, &
                          test3allB => eos_test3allB, &
                          test4allB => eos_test4allB
  implicit none

# include "Eos.h"
# include "constants.h"
# include "Flash.h"

  integer, intent(in) :: fileUnit
  logical, intent(out) :: perfect
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t),OPTIONAL, intent(in) :: blockDesc

  integer :: blockID
  real :: presErr, tempErr, eintErr

  logical:: test1,test2,test3,test4 !for a block

  integer :: vecLen, blockOffset,  pres, dens, temp, e, n, m
  integer :: isize, jsize, ksize, i,j,k, nStartsAtOne
  real, dimension(:), allocatable :: eosData
  real, dimension(:), allocatable :: massFrac
  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
  real, allocatable, dimension(:,:,:,:) :: derivedVariables

  character(len=7),pointer:: ap
  character(len=7),target :: a
  integer,parameter :: maxPrintPE = 20
  integer :: level
  integer :: nodeType
  integer :: ib,ie,jb,je,kb,ke
  real presErr1, presErr2


  if (eos_meshNumProcs==1) then
     a = ''
     ap => a
  else
20   format(I6,':')
     write(a,20) eos_meshMe
     a = trim(adjustl(a))
     ap => a
  end if

! info for checkpoint output


  mask = .true.

  !! In Simulation_initBlock,
  !! temperature is initialized in CTMP_VAR and pressure is
  !! initialized in CPRS_VAR. We don't change these variables
  !! we copy them into the usual variable name as needed.
  
  ib=blkLimits(LOW,IAXIS)
  ie=blkLimits(HIGH,IAXIS)
  
  jb=blkLimits(LOW,JAXIS)
  je=blkLimits(HIGH,JAXIS)
  
  kb=blkLimits(LOW,KAXIS)
  ke=blkLimits(HIGH,KAXIS)

  ! This is just to get more information into the output, FLASH4-like
  blockID = -1; nodeType = -1; level = -1
  if (present(blockDesc)) then
#ifdef FLASH_GRID_PARAMESH
     blockID = blockDesc%id
     call Grid_getBlkType(blockID, nodeType)
#else
     blockID = blockDesc % grid_index
     nodeType = -1                !We do not really know.
#endif
     level = blockDesc % level
  end if

  ! Testing density/temperature in; energy/pressure out
  if (eos_meshMe<maxPrintPE .AND. level > -1) print *,ap,'Block',blockID,' type',nodeType,' level',level
  if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_TEMP or similar: Density, temperature in; energy, pressure out; mode=', &
       eos_testTempMode,eos_testTempModeStr
  if (eos_meshMe<maxPrintPE) then
     print*,ap,'The initialized extreme values are '
     print*,ap,'Initialized Density min abs',minval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Initialized Density max abs',maxval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Initialized Temperature min abs',minval(abs(solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Initialized Temperature max abs',maxval(abs(solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Initialized Pressure min',minval(solnData(CPRS_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Initialized Pressure max',maxval(solnData(CPRS_VAR,ib:ie,jb:je,kb:ke))
  end if
  
  solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)
  call Eos_wrapped(eos_testTempMode, blkLimits,solnData, CENTER)
  !! Summarize results of MODE_DENS_TEMP (or similar) call
  if (eos_meshMe<maxPrintPE) then
     print*,ap,'The resulting extreme values are '
        print*,ap,'Resulting Density min',minval(solnData(DENS_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'Resulting Density max',maxval(solnData(DENS_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting Temperature min',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting Temperature max',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting Pressure min',minval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting Pressure max',maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting internal energy min',minval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting internal energy max',maxval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
  end if
  
  ! Save the equilibrium values in O variables and (hopefully) write them out
  solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
  solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
  solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)
  


!!$  test1allB = .TRUE.
  !  Testing density/energy in, temperature/pressure out
  if (eos_meshMe<maxPrintPE .AND. level > -1) print *,ap,'Block',blockID,' type',nodeType,' level',level
  if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI or similar: Density, energy in; temperature, pressure out; mode=', &
       eos_testEintMode,eos_testEintModeStr
  !  Zero output variables
  solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=1.e-10  ! don't zero TEMP or eos_helm cannot converge in MODE_DENS_EI
  solnData(PRES_VAR,:,:,:)=0 
  call Eos_wrapped(eos_testEintMode,blkLimits,solnData,CENTER)

  
  if (eos_meshMe<maxPrintPE) then !! Summarize results of MODE_DENS_EI (or similar) call
     print*,ap,'  Temperature min ',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'  Temperature max ',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'  Pressure min ',minval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'  Pressure max ',maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
  end if
  
  !! Calculate error from MODE_DENS_EI (or similar) call.
  tempErr = maxval(   &
       abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)- solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/ &
       solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))  )
  if (eos_meshMe<maxPrintPE) then
     print*,ap,'  The calculated error in temperature is ',tempErr
  end if
  
  presErr = maxval(abs((solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/&
       solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
  if (eos_meshMe<maxPrintPE) print*,ap,'  The calculated error in pressure is ',presErr
  

  test1 = (tolerance > tempErr)
  test1 = test1.and.(tolerance > presErr)
  if(test1) then
     if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_EI or similar is fine'
  else
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI or similar is BAD!!!'
     test1allB = .FALSE.
  endif




!!$  test2allB = .TRUE.
! Testing density/pressure in, energy/temp out
  if (eos_meshMe<maxPrintPE .AND. level > -1) print *,ap,'Block',blockID,' type',nodeType,' level',level
  if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_PRES or similar: Density, pressure in; energy, temperature out; mode=', &
       eos_testPresMode,eos_testPresModeStr
  solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0
  solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=1.1e4  ! don't zero TEMP or eos_helm cannot converge in any mode
  call Eos_wrapped(eos_testPresMode,blkLimits,solnData,CENTER)

  !! Summarize results of MODE_DENS_PRES (or similar) call;
  !! calculate error from MODE_DENS_PRES (or similar) call.
  tempErr = maxval(abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)-&
       solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
  eintErr = maxval(abs((solnData(EINT_VAR,ib:ie,jb:je,kb:ke)-&
       solnData(OENT_VAR,ib:ie,jb:je,kb:ke))/solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
  if (eos_meshMe<maxPrintPE) print*,ap,'  Energy min is',minval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
  if (eos_meshMe<maxPrintPE) print*,ap,'  Energy max is',maxval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
  if (eos_meshMe<maxPrintPE) print*,ap,'  Temperature min is',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
  if (eos_meshMe<maxPrintPE) print*,ap,'  Temperature max is',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
  if (eos_meshMe<maxPrintPE) print*,ap,'  The calculated error in energy is ',eintErr
  if (eos_meshMe<maxPrintPE) print*,ap,'  The calculated error in temperature is ',tempErr
  
  test2 = (tolerance > tempErr)
  test2 = test2.and.(tolerance > eintErr)
  if(test2) then
     if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_PRES or similar is fine'
  else
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_PRES or similar is BAD!!!'
     test2allB = .FALSE.
  endif
  

!!$  test3allB = .TRUE.
!!$  test4allB = .TRUE.
  if (eos_meshMe<maxPrintPE .AND. level > -1) print *,ap,'Block',blockID,' type',nodeType,' level',level
  if (eos_meshMe<maxPrintPE) print*,ap,'And now to verify the other solutions'
  
  ! Generate the initial conditions again
  solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=solnData(CPRS_VAR,ib:ie,jb:je,kb:ke)
  solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)
  
  !!  Do calculations in reverse order
  ! Density and pressure in, energy and temperature out
  !solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=0   ! don't zero TEMP or eos_helm cannot converge
  solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0 
  call Eos_wrapped(MODE_DENS_PRES, blkLimits,solnData)
  if (eos_meshMe<maxPrintPE) then
     print*,ap,'The resulting extreme values from MODE_DENS_PRES are '
     print*,ap,'Resulting Pressure min',minval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting Pressure max',maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting Temperature min abs',minval(abs(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Resulting Temperature max abs',maxval(abs(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Resulting internal energy min abs',minval(abs(solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Resulting internal energy max abs',maxval(abs(solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'Resulting Temperature min',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting Temperature max',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting internal energy min',minval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
     print*,ap,'Resulting internal energy max',maxval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
     presErr = maxval(abs((solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-solnData(CPRS_VAR,ib:ie,jb:je,kb:ke))/&
          solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
     print*,ap,'The calculated error in pressure is ',presErr
  end if
  ! Now we have a "true"  temperature and internal energy; save them for comparison
  solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
  solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
  solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)
  
  
  ! Density and energy in, temperature and pressure out
  !! zero output values to make sure they're being calculated
  solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=0.0
  !solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=1.0e-11   ! don't zero TEMP or eos_helm cannot converge
  call Eos_wrapped(MODE_DENS_EI,blkLimits,solnData,CENTER)
  presErr1 = maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
  presErr2 = maxval(solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))
  if (eos_meshMe<maxPrintPE) print *,ap,'maxval PRES_VAR OPRS_VAR',presErr1,presErr2
  presErr = maxval(abs(&
       (solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/&
       solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
  ! NOTE  this ALWAYS comes out to zero. suspect  something wrong here....
  tempErr = maxval(abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)-&
       solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
  test3 = (tolerance > tempErr)
  test3 = test3.and.(tolerance > presErr)
  if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in pressure from EI is ',presErr
  if (eos_meshMe<maxPrintPE) print*,ap,'  Temperature min is',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
  if (eos_meshMe<maxPrintPE) print*,ap,'  Temperature max is',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
  if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in temperature from EI is ',tempErr
  if(test3) then
     if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_EI is fine '
  else
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI is BAD!!!'
     test3allB = .FALSE.
  endif
  
  solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0
  solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=0 
  call Eos_wrapped(MODE_DENS_TEMP,blkLimits,solnData,CENTER)
  presErr = maxval(abs((solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-&
       solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
  eintErr = maxval(abs((solnData(EINT_VAR,ib:ie,jb:je,kb:ke)-&
       solnData(OENT_VAR,ib:ie,jb:je,kb:ke))/solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
  test4 = (tolerance > presErr)
  test4 = test4.and.(tolerance > eintErr)
  if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in pressure is ',presErr
  if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in energy is ',eintErr
  if(test4) then
     if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_TEMP is fine'
  else
     if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_TEMP is BAD!!!!'
     test4allB = .FALSE.
  endif
  solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
  solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
  solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)
  
  !! Finally, do a test of the derived variables just for exercise.....
  if (eos_meshMe<maxPrintPE) print *,ap,' Now testing the derived variables'
  
  !  Allocate the necessary arrays for an entire block of data
  isize = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
  jsize = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
  ksize = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)
  
  vecLen=isize
  
  allocate(derivedVariables(isize,jsize,ksize,EOS_NUM))
  allocate(eosData(vecLen*EOS_NUM))
  allocate(massFrac(vecLen*NSPECIES))
  
  ! Initialize them
  derivedVariables = 0.0
  vecLen=isize
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  
     
  !! Get DENS and PRES to fill up input, also massFraction
  solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
        do i = 1,vecLen
           massFrac((i-1)*NSPECIES+1:i*NSPECIES) = &
                solnData(SPECIES_BEGIN:SPECIES_END,ib+i-1,j,k)
        end do
        
        eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
        eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
        eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
        
        call Eos(MODE_DENS_PRES,vecLen,eosData,massFrac,mask)
           
        do e=EOS_VARS+1,EOS_NUM
           m = (e-1)*vecLen
           derivedVariables(1:vecLen,j-jb+1,k-kb+1,e) =  eosData(m+1:m+vecLen)
        end do
     end do
  end do
  
  
  
  if (All(derivedVariables .eq. 0.0)) then
     if (eos_meshMe<maxPrintPE) print*,ap,"No derived variables were set!"
  end if
  if (ANY(derivedVariables .ne. 0.0)) then
     if (eos_meshMe<maxPrintPE) print*,ap,"Some derived variables were set."
  end if
  deallocate(eosData)
  deallocate(massFrac)
  deallocate(derivedVariables)
  
  perfect = test1allB.and.test2allB.and.test3allB.and.test4allB
  if(perfect) then
     if (eos_meshMe<maxPrintPE) print*,ap,'SUCCESS(so far) all tests were fine'
  else
     if (eos_meshMe<maxPrintPE) then
        print*,ap,'FAILURE(so far) some tests failed'
        if (.NOT.(test1.and.test2.and.test3.and.test4)) then
           print*,ap,'FAILURE(block) tests1..4:',test1,test2,test3,test4
        end if
     end if
  end if
  return
end subroutine Eos_unitTest4




