/*@@
  @file      setupbrilldata3D.F
  @date      October 1998
  @author    Miguel Alcubierre
  @desc
             Set up non-axisymmetric Brill data for elliptic solve.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine setupbrilldata3D(CCTK_ARGUMENTS)

c     Set up 3D Brill data for elliptic solve.  The elliptic
c     equation we need to solve is:
c
c     __
c     \/  psi  -  psi R  /  8  =  0
c       c              c
c
c     where:
c     __
c     \/  =  Laplacian operator for conformal metric.
c       c
c
c     R   =  Ricci scalar for conformal metric. 
c      c
c
c     The Ricci scalar for the conformal metric turns out to be:
c
c               /  -2q    2      2             -2            2       2      \
c     R  =  - 2 | e    ( d q  + d   q )  +  rho   ( 3 (d   q)  +  2 d   q ) |
c      c        \         z      rho                    phi          phi    /


      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer CCTK_Equals
      integer i,j,k
      integer nx,ny,nz
      integer ierr

      integer, dimension(3) :: sym

      CCTK_REAL x1,y1,z1,rho1,rho2
      CCTK_REAL phi,phif,e2q
      CCTK_REAL brillq,rhofudge,eps
      CCTK_REAL zero,one,two,three

c     Set up grid size.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

c     Define numbers

      zero  = 0.0D0
      one   = 1.0D0
      two   = 2.0D0
      three = 3.0D0

c     Parameters.

      rhofudge = brill_rhofudge

c     Epsilon for finite differencing.

      eps = cctk_delta_space(1)

c     Initialize psi.

      brillpsi = one

c     Set up conformal metric.

      if (use_conformal == 1) then
      
         conformal_state = CONFORMAL_METRIC
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  psi(i,j,k) = one
                  psix(i,j,k) = zero
                  psiy(i,j,k) = zero
                  psiz(i,j,k) = zero
                  psixx(i,j,k) = zero
                  psiyy(i,j,k) = zero
                  psizz(i,j,k) = zero
                  psixy(i,j,k) = zero
                  psixz(i,j,k) = zero
                  psiyz(i,j,k) = zero
               end do
            end do
         end do

      else

         conformal_state = NOCONFORMAL_METRIC

      end if
         

      do k=1,nz
         do j=1,ny
            do i=1,nx

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho2 = x1**2 + y1**2
               rho1 = dsqrt(rho2)

               phi = phif(x1,y1)

               e2q  = dexp(two*brillq(rho1,z1,phi))

               if (rho1.gt.rhofudge) then
                  gxx(i,j,k) = e2q + (one - e2q)*y1**2/rho2
                  gyy(i,j,k) = e2q + (one - e2q)*x1**2/rho2
                  gzz(i,j,k) = e2q
                  gxy(i,j,k) = - (one - e2q)*x1*y1/rho2
               else
                  gxx(i,j,k) = one 
                  gyy(i,j,k) = one 
                  gzz(i,j,k) = one
                  gxy(i,j,k) = zero
               end if

            end do
         end do
      end do

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

               rho2 = x1**2 + y1**2
               rho1 = dsqrt(rho2)

               phi  = phif(x1,y1)

               e2q = dexp(two*brillq(rho1,z1,phi))

c              Find M using centered differences over a small
c              interval.

               if (rho1.gt.rhofudge) then

                  brillMlinear(i,j,k) = 0.25D0/e2q
     .                 *(brillq(rho1,z1+eps,phi) 
     .                 + brillq(rho1,z1-eps,phi) 
     .                 + brillq(rho1+eps,z1,phi) 
     .                 + brillq(rho1-eps,z1,phi) 
     .                 - 4.0D0*brillq(rho1,z1,phi))
     .                 / eps**2

               else

                  brillMlinear(i,j,k) = 0.25D0/e2q
     .                 *(brillq(rho1,z1+eps,phi) 
     .                 + brillq(rho1,z1-eps,phi) 
     .                 + two*brillq(eps,z1,phi) 
     .                 - two*brillq(rho1,z1,phi))
     .                 / eps**2

               end if

c              Here I assume that for very small rho, the
c              phi derivatives are essentially zero.  This
c              must always be true otherwise the function
c              is not regular on the axis.

               if (rho1.gt.rhofudge) then
                  brillMlinear(i,j,k) = brillMlinear(i,j,k) + 0.25D0/rho2
     .               *(three*0.25D0*(brillq(rho1,z1,phi+eps)
     .               - brillq(rho1,z1,phi-eps))**2
     .               + two*(brillq(rho1,z1,phi+eps)
     .               - two*brillq(rho1,z1,phi)
     .               + brillq(rho1,z1,phi-eps)))/eps**2
               end if

            end do
         end do
      end do

c     Set coefficient Nsource = 0.

      brillNsource = zero

c     Synchronization and boundaries.

      call CCTK_SyncGroup(ierr,cctkGH,'idbrilldata::brillelliptic')
      call CartSymGN(ierr,cctkGH,'idbrilldata::brillelliptic')

      return
      end

