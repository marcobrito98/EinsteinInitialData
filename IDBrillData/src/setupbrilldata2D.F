/*@@
  @file      setupbrilldata2D.F
  @date
  @author    Carsten Gundlach (Cactus 4, Miguel Alcubierre)
  @desc
             Set up axisymmetric Brill data for elliptic solve.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine setupbrilldata2D(CCTK_FARGUMENTS)

c     Set up axisymmetric Brill data for elliptic solve. The elliptic
c     equation we need to solve is:
c
c     __                 2      2
c     \/  psi  +  psi ( d q  + d   q ) / 4  =  0
c       f                z      rho
c
c     where:
c
c     __
c     \/  =  Flat space Laplacian.
c       f

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k
      integer nx,ny,nz
      integer ierr

      integer, dimension(3) :: sym

      CCTK_REAL x1,y1,z1,rho1
      CCTK_REAL brillq,eps
      CCTK_REAL zero,one

      external brillq

c     Numbers.

      zero = 0.0D0
      one  = 1.0D0

c     Set up grid size.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

c     Parameters.

      eps = brill_eps

c     Initialize psi.

      brillpsi = one

c     Initialize metric.

      psi = one

      gxx = one
      gyy = one
      gzz = one

      gxy = zero
      gxz = zero
      gyz = zero

c     Define the symmetries for the brill GFs.

      sym(1) = 1
      sym(2) = 1
      sym(3) = 1

      call SetCartSymVN(ierr,cctkGH,sym,'idbrilldata::brillMlinear')
      call SetCartSymVN(ierr,cctkGH,sym,'idbrilldata::brillNsource')

c     Set up coefficient arrays for elliptic solve.
c     Notice that the Cactus conventions are:
c     __
c     \/ psi +  Mlinear*psi  +  Nsource  =  0

      do k=1,nz
         do j=1,ny
            do i=1,nx

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho1 = dsqrt(x1**2 + y1**2)

c              Centered derivatives. Note that here we may be calling brillq
c              with a small negative rho, but that should be ok as long as
c              brillq is even in rho - physically it must be, or the data
c              will not be regular on the axis.

               brillMlinear(i,j,k) = 0.25d0 
     .              *(brillq(rho1,z1+eps,zero) 
     .              + brillq(rho1,z1-eps,zero) 
     .              + brillq(rho1+eps,z1,zero) 
     .              + brillq(rho1-eps,z1,zero) 
     .              - 4.d0*brillq(rho1,z1,zero))/eps**2

            end do
         end do
      end do

c     Set coefficient Nsource = 0.

      brillNsource = zero

c     Synchronization and boundaries.

      call CCTK_SyncGroup(cctkGH,'idbrilldata::brillelliptic')
      call CartSymGN(ierr,cctkGH,'idbrilldata::brillelliptic')

      return
      end
