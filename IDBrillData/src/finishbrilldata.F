#include "cctk.h" 
#include "cctk_parameters.h"
#include "cctk_arguments.h"
#include "../../../CactusEinstein/Einstein/src/Einstein.h"

      subroutine finishbrilldata(CCTK_FARGUMENTS)

C     Author: Carsten Gundlach.
C
C     Set up 3-metric, extrinsic curvature and BM variables
C     once the conformal factor has been found.

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k
      integer nx,ny,nz
      integer CCTK_Equals

      CCTK_REAL x1,y1,z1,rho1,rho2
      CCTK_REAL brillq,psi4,e2q,rhofudge
      CCTK_REAL phi,phif
      CCTK_REAL zero

      external brillq

C     Numbers.

      zero = 0.0D0

C     Set up grid size.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     Parameters.

      rhofudge = brill_rhofudge
      
C     Replace flat metric left over from elliptic solve by
C     Brill metric calculated from q and Psi.

      do k=1,nz
         do j=1,ny
            do i=1,nx

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho2 = x1**2 + y1**2
               rho1 = dsqrt(rho2)

               phi = phif(x1,y1)

               psi4 = brillpsi(i,j,k)**4
               e2q  = dexp(2.d0*brillq(rho1,z1,phi))

C              Fudge division by rho^2 on axis. (Physically, y^/rho^2,
C              x^2/rho^2 and xy/rho^2 are of course regular.)
C              Transform Brills form of the physical metric to Cartesian
C              coordinates via
C
C              e^2q (drho^2 + dz^2) + rho^2 dphi^2 =
C              e^2q (dx^2 + dy^2 + dz^2) + (1-e^2q) (xdy-ydx)^2/rho^2
C
C               The individual coefficients can be read off as

               if (rho1 .gt. rhofudge) then

                  gzz(i,j,k) = psi4*e2q
                  gxx(i,j,k) = psi4*(e2q + (1.d0 - e2q)*y1**2/rho2)
                  gyy(i,j,k) = psi4*(e2q + (1.d0 - e2q)*x1**2/rho2)
                  gxy(i,j,k) = - psi4*(1.d0 - e2q)*x1*y1/rho2

               else

C              This fudge assumes that q = O(rho^2) near the axis. Which
C              it should be, or the data will be singular.

                  gzz(i,j,k) = psi4 
                  gxx(i,j,k) = psi4 
                  gyy(i,j,k) = psi4 
                  gxy(i,j,k) = zero

               end if

            end do
         end do
      end do

C     In any case,

      gxz = zero
      gyz = zero

C     Vanishing extrinsic curvature completes the Cauchy data.

      kxx = zero
      kyy = zero
      kzz = zero
      kxy = zero
      kxz = zero
      kyz = zero

      return
      end


