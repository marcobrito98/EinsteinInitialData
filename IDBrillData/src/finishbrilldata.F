/*@@
  @file      finishbrilldata.F
  @date
  @author    Carsten Gundlach (Cactus 4, Miguel Alcubierre)
  @desc
             Reconstruct metric from conformal factor.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine finishbrilldata(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k
      integer nx,ny,nz

      CCTK_REAL x1,y1,z1,rho1,rho2
      CCTK_REAL phi,psi4,e2q,rhofudge
      CCTK_REAL zero,one

      CCTK_REAL brillq,phif

c     Numbers.

      zero = 0.0D0
      one  = 1.0D0

c     Set up grid size.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

c     Parameters.

      rhofudge = brill_rhofudge
      
c     Replace flat metric left over from elliptic solve by
c     Brill metric calculated from q and Psi.

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

c              Fudge division by rho^2 on axis. (Physically, y^/rho^2,
c              x^2/rho^2 and xy/rho^2 are of course regular.)
c              Transform Brills form of the physical metric to Cartesian
c              coordinates via
c
c              e^2q (drho^2 + dz^2) + rho^2 dphi^2 =
c              e^2q (dx^2 + dy^2 + dz^2) + (1-e^2q) (xdy-ydx)^2/rho^2
c
c              The individual coefficients can be read off as

               if (rho1.gt.rhofudge) then

                  gxx(i,j,k) = psi4*(e2q + (one - e2q)*y1**2/rho2)
                  gyy(i,j,k) = psi4*(e2q + (one - e2q)*x1**2/rho2)
                  gzz(i,j,k) = psi4*e2q
                  gxy(i,j,k) = - psi4*(one - e2q)*x1*y1/rho2

               else

c                 This fudge assumes that q = O(rho^2) near the axis. Which
c                 it should be, or the data will be singular.

                  gxx(i,j,k) = psi4
                  gyy(i,j,k) = psi4
                  gzz(i,j,k) = psi4
                  gxy(i,j,k) = zero

               end if

            end do
         end do
      end do

c     In any case,

      gxz = zero
      gyz = zero

c     Vanishing extrinsic curvature completes the Cauchy data.

      kxx = zero
      kyy = zero
      kzz = zero
      kxy = zero
      kxz = zero
      kyz = zero

      return
      end


