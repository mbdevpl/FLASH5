!!****if* source/Simulation/SimulationComposition/Burn/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed
!!  for setups that use nuclear networks.  The setups that want to use multispecies
!!  capabilities of the code for something other than nuclear burning and ionization
!!  should include their own custom implementation of this routine
!!
!!  This implementation of the routine relies on a textfile
!!  SpeciesList.txt to provide the species related informatio. The
!!  textfile contains the elements sorted by their atomic number in
!!  increasing order, and the isotopes of each element in turn sorted by
!!  of their atomic number, again in increasing order. The subroutine
!!  reads in the records in file, if the record corresponds to an
!!  isotope that is included in the setup, it sets the properties of
!!  the isotope in the multispecies database, and if the istope is not
!!  included, it goes on to read the next one. This process is
!!  repeated until all the species included in the setup have been
!!  found.
!!
!!  The format of SpeciesList.txt is as follow
!!  Column#         Variable                  Description
!!  ---------------------------------------------------------
!!  1               isotopeName               Sorted in increasing atomic number
!!  2               Z                         zbar, Atomic number; number of protons in nucleus
!!  3               A                         abar, total number of protons and neutrons in nucleus
!!  4               N                         Number of neutrons
!!  5               Eb                        binding energy
!!  6               ??
!!  7               ??
!!
!!
!!
!!  ARGUMENTS : There are no arguments in this subroutine
!!
!!  NOTE
!! 
!!***

subroutine Simulation_initSpecies()
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_setProperty
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use xnet_constants, ONLY : m_e, m_n, m_p, m_u
  implicit none

#include "constants.h"
#include "Flash.h"
#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
  
  character(len=80) :: data_dir
  character(len=4) :: isotopeName
  character(len=5) :: nname(NSPECIES), tmp_name
  real, dimension(NSPECIES) :: aa, zz, nn, be, mex
  real :: spin, mex_n, mex_p
  integer :: ny, iz, in, ineut, iprot, lun_winv
  integer :: i, j, isotope

  ! Get the species data from XNet data files
  call RuntimeParameters_get('xnet_data_dir',data_dir)
  open(newunit=lun_winv, file=trim(data_dir)//"/netwinv", status='old')
  read(lun_winv,"(i5)") ny
  if ( ny /= NSPECIES ) then
     call Driver_abortFlash("[Simulation_initSpecies] ny /= NSPECIES")
  end if

  ! Skip the header information
  read(lun_winv,*)
  do i = 1, ny
     read(lun_winv,*)
  end do
  
  ineut = 0 ; iprot = 0
  do i = 1, ny
     read(lun_winv,*) nname(i), aa(i), iz, in, spin, mex(i)
     do j = 1, 3
        read(lun_winv,*)
     end do
     if ( iz == 0 .and. in == 1 ) ineut = i
     if ( iz == 1 .and. in == 0 ) iprot = i
     zz(i) = real(iz)
     nn(i) = real(in)
  end do
  close(lun_winv)

  ! For consistency, use neutron and proton mass excess from netwinv to
  ! calculate binding energies. Otherwise, use CODATA recommended 2014 values.
  if ( ineut > 0 ) then
     mex_n = mex(ineut)
  else
     mex_n = m_n - m_u
  end if
  if ( iprot > 0 ) then
     mex_p = mex(iprot)
  else
     mex_p = m_p + m_e - m_u
  end if
  be(:) = mex_n*nn(:) + mex_p*zz(:) - mex(:)

  do i = 1, NSPECIES
     tmp_name = adjustl(nname(i))
     isotopeName = tmp_name(1:4)
     call Simulation_mapStrToInt(isotopeName,isotope,MAPBLOCK_UNK)
     if ( isotope /= NONEXISTENT ) then
        call Multispecies_setProperty(isotope, A, aa(i))
        call Multispecies_setProperty(isotope, Z, zz(i))
        call Multispecies_setProperty(isotope, N, nn(i))
        call Multispecies_setProperty(isotope, E, zz(i))
        call Multispecies_setProperty(isotope, EB, be(i))
     else
        call Driver_abortFlash("[Simulation_initSpecies] Isotope missing from UNK: "//isotopeName)
     end if
  end do
#endif

end subroutine Simulation_initSpecies
