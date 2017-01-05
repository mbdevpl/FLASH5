!!****ih* source/diagnostics/ThomsonScattering/localAPI/thsc_constInterface
!!
!! NAME
!!
!!  thsc_constInterface
!!
!! SYNOPSIS
!!
!!   use thsc_constInterface
!!
!! DESCRIPTION
!!
!!  This is a private header file for the Thomson Scattering unit that defines
!!  some constants.
!!
!!***

Module thsc_constInterface

  ! for integration results: (optical depth, Faraday rotation angle, etc.)
  integer,parameter :: RAYSUM_SIZE    = 7, &
                       RAYSUM_CRITIND = 1, &
                       RAYSUM_DEPTH   = 2, &
                       RAYSUM_FARRO   = 3, &
                       RAYSUM_TIME    = 4, &
                       RAYSUM_ATIME   = 5, &
                       RAYSUM_ALEN    = 6, &
                       RAYSUM_LEN     = 7

end Module thsc_constInterface
