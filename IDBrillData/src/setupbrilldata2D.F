
#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine setupbrilldata2D(CCTK_FARGUMENTS)

C     Author: Carsten Gundlach.
C
C     Set up axisymmetric Brill data for elliptic solve. The elliptic
C     equation we need to solve is:
C
C     __                 2      2
C     \/  psi  =  psi ( d q  + d   q ) / 4
C       f                z      rho
C
C     where:
C
C     __
C     \/  =  Flat space Laplacian.
C       f

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k
      integer nx,ny,nz

      CCTK_REAL x1,y1,z1,rho1
      CCTK_REAL brillq,eps
      CCTK_REAL zero,one

      external brillq

C     Numbers.

      zero = 0.0D0
      one  = 1.0D0

C     Set up grid size.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     Parameters.

      eps = brill_eps

C     Initialize psi.

      brillpsi = one

C     Initialize metric.

      psi = one

      gxx = one
      gyy = one
      gzz = one

      gxy = zero
      gxz = zero
      gyz = zero

C     Set up coefficient Mlinear = - (1/4) (q,zz + q,rhorho).

      do k=1,nz
         do j=1,ny
            do i=1,nx

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho1 = dsqrt(x1**2 + y1**2)

C              Centered derivatives. Note that here we may be calling brillq
C              with a small negative rho, but that should be ok as long as
C              brillq is even in rho - physically it must be, or the data
C              will not be regular on the axis.

               brillMlinear(i,j,k) = - 0.25d0 
     $              *(brillq(rho1,z1+eps,zero) 
     $              + brillq(rho1,z1-eps,zero) 
     $              + brillq(rho1+eps,z1,zero) 
     $              + brillq(rho1-eps,z1,zero) 
     $              - 4.d0*brillq(rho1,z1,zero))/eps**2

            end do
         end do
      end do

C     Set coefficient Nsource = 0.

      brillNsource = zero

      return
      end
