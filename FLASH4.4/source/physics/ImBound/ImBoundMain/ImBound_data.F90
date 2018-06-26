!!****if* source/physics/ImBound/ImBoundMain/ImBound_data
!!
!! NAME
!!
!!  ImBound_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE ImBound_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data specific to the Immersed Boundary module.
!!
!!***


module ImBound_data
     
      real, save :: ib_nu = 0. 
      real, save :: ib_dt      

      integer,save :: ib_meshMe  
      integer,save :: ib_globalMe 

end module ImBound_data

