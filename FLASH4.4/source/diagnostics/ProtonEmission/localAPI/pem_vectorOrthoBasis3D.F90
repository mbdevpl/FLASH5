!!****if* source/diagnostics/ProtonEmission/localAPI/pem_vectorOrthoBasis3D
!!
!! NAME
!!
!!  pem_vectorOrthoBasis3D
!!
!! SYNOPSIS
!!
!!  call pem_vectorOrthoBasis3D (real, intent (in)  :: Vx,
!!                               real, intent (in)  :: Vy,
!!                               real, intent (in)  :: Vz,
!!                               real, intent (out) :: ix,
!!                               real, intent (out) :: iy,
!!                               real, intent (out) :: iz,
!!                               real, intent (out) :: jx,
!!                               real, intent (out) :: jy,
!!                               real, intent (out) :: jz,
!!                               real, intent (out) :: kx,
!!                               real, intent (out) :: ky,
!!                               real, intent (out) :: kz)
!!
!! DESCRIPTION
!!
!!  Given a vector V with 3D cartesian components (Vx,Vy,Vz), this routine establishes
!!  an orthogonal basis i,j,k at the vector origin, such that k is collinear with V
!!  and i,j,k form a right handed system. Hence the following relations are fulfilled:
!!
!!                               |i| , |j| , |k| = 1     magnitudes
!!                               i.j , i.k , j.k = 0     scalar products
!!                                         i x j = k     cross product
!!                                     i.V , j.V = 0     scalar products
!!                                           k.V = |V|   scalar product
!!
!!  Although there are an infinite set of orthogonal bases due to possible rotation
!!  of the i,j pair around the V axis, we choose the one that has its most simple form,
!!  of which there are three. The first such case is:
!!
!!                            ix = + |Vyz| / |V|
!!                            iy = - Vx * Vy / (|Vyz| * |V|)
!!                            iz = - Vx * Vz / (|Vyz| * |V|)
!!                            jx = 0
!!              Case I        jy = + Vz / |Vyz|
!!                            jz = - Vy / |Vyz|
!!                            kx = Vx / |V|
!!                            ky = Vy / |V|
!!                            kz = Vz / |V|
!!
!!  where |V| = sqrt (Vx*Vx + Vy*Vy + Vz*Vz) and |Vyz| = sqrt (Vy*Vy + Vz*Vz). Computationally
!!  however, it is best to calculate these as the following sequence:
!!
!!                            kx = Vx / |V|
!!                            ky = Vy / |V|
!!                            kz = Vz / |V|
!!                            ix = sqrt (ky * ky + kz * kz)
!!              Case I        jx = 0
!!                            jy = kz / ix
!!                            jz = - ky / ix
!!                            iy = jz * kx
!!                            iz = - jy * kx
!!
!!  The other remaining 2 cases are obtained by performing the label permutations:
!!
!!              Case II       x -> y , y -> z , z -> x
!!              Case III      x -> z , y -> x , z -> y
!!
!!  We choose among the cases as follows, depending on the magnitude of the kx,ky,kz:
!!
!!                         |ky|,|kz| largest -> Case I
!!                         |kx|,|kz| largest -> Case II
!!                         |kx|,|ky| largest -> Case III
!!
!! ARGUMENTS
!!
!!  Vx : x-coordinate of vector V
!!  Vy : y-coordinate of vector V
!!  Vz : z-coordinate of vector V
!!  ix : x-coordinate of orthogonal basis vector i
!!  iy : y-coordinate of orthogonal basis vector i
!!  iz : z-coordinate of orthogonal basis vector i
!!  jx : x-coordinate of orthogonal basis vector j
!!  jy : y-coordinate of orthogonal basis vector j
!!  jz : z-coordinate of orthogonal basis vector j
!!  kx : x-coordinate of orthogonal basis vector k
!!  ky : y-coordinate of orthogonal basis vector k
!!  kz : z-coordinate of orthogonal basis vector k
!!
!! NOTES
!!
!!  Orthogonality and normality is not tested.
!!
!!***

subroutine pem_vectorOrthoBasis3D (Vx,Vy,Vz,  ix,iy,iz,jx,jy,jz,kx,ky,kz)

  implicit none

  real, intent (in)  :: Vx, Vy, Vz
  real, intent (out) :: ix, iy, iz
  real, intent (out) :: jx, jy, jz
  real, intent (out) :: kx, ky, kz

  ix = 0.0
  iy = 0.0
  iz = 0.0
  jx = 0.0
  jy = 0.0
  jz = 0.0
  kx = 0.0
  ky = 0.0
  kz = 0.0

  return
end subroutine pem_vectorOrthoBasis3D
