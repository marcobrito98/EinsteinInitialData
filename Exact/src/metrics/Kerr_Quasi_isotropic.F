c Kerr metric in Quasi-isotropic coordinates,
c as per arXiv:1410.0694v2
C
C Author: Marco Brito
C Copyright/License: unknown
C
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Kerr_Quasi_isotropic(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz,
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz,
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

         implicit none

         DECLARE_CCTK_PARAMETERS

c input arguments
         CCTK_REAL x, y, z
         CCTK_DECLARE(CCTK_REAL, t,)

c output arguments
         CCTK_REAL gdtt, gdtx, gdty, gdtz,
     $          gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $          gutt, gutx, guty, gutz,
     $          guxx, guyy, guzz, guxy, guyz, guzx
         CCTK_DECLARE(CCTK_REAL, psi,)
         LOGICAL   Tmunu_flag

c local static variables
         CCTK_REAL eps, bhmass, bhspin, posx, posy, posz

c local variables
         CCTK_REAL r, r2, rho2, costh, costh2, sinth, sinth2
         CCTK_REAL horizonradius, rBL, rho2kerr, sigma, hh, psi4
         CCTK_REAL alpha0, bphi, gammaphiphi, bphiup, betadX, betadY
         CCTK_REAL x1, y1, z1
         CCTK_REAL gdn(0:3,0:3), gun(0:3,0:3), prod(0:3,0:3)
         CCTK_REAL tol, maxerr, err
         INTEGER i, j, k
         CHARACTER(len=256) msg

C This is a vacuum spacetime with no cosmological constant
         Tmunu_flag = .false.

C     Get parameters of the exact solution.

         eps = Kerr_KerrQI__eps_r
         bhmass = Kerr_KerrQI__mass
         bhspin = Kerr_KerrQI__spin
         posx = Kerr_KerrQI__x
         posy = Kerr_KerrQI__y
         posz = Kerr_KerrQI__z

         x1 = x - posx
         y1 = y - posy
         z1 = z - posz

         r2 = max(x1**2 + y1**2 + z1**2, eps**2)
         r = sqrt(r2)
         rho2 = max(x1**2 + y1**2, eps**2)
         costh = z1 / r
         costh2 = costh**2
         sinth2 = 1.d0 - costh2
         sinth = sqrt(sinth2)

         horizonradius = 0.5 * sqrt(bhmass**2 - bhspin**2)
         rBL = r + bhmass + 0.25*(bhmass**2 - bhspin**2) / r
         rho2kerr = rBL**2 + bhspin**2 * costh2
         sigma = 2.d0 * bhmass * rBL / rho2kerr
         hh = (1.d0 + sigma) / (rho2kerr * r2)
         psi4 = rho2kerr / r2
         alpha0 = (r + horizonradius) * (r - horizonradius) / r *
     $   1.d0 / sqrt(rBL**2 + bhspin**2 * (1.d0 + sigma * sinth2))
         bphi = - bhspin * sigma * sinth2
         gammaphiphi = rho2kerr * rho2 / r2 * (1.d0 +
     $    bhspin**2 * hh * rho2)
         bphiup = bphi / gammaphiphi
         betadX = -y1 /rho2 * bphi
         betadY =  x1 /rho2 * bphi

         gdtt = - alpha0**2 + bphi*bphiup
         gdtx = betadX
         gdty = betadY
         gdtz = 0.d0
         gdxx = psi4 * (1.d0 + bhspin**2 * hh * y1**2)
         gdyy = psi4 * (1.d0 + bhspin**2 * hh * x1**2)
         gdzz = psi4
         gdxy = - psi4 * bhspin**2 * hh * y1 * x1
         gdyz = 0.d0
         gdzx = 0.d0

         gutt = -(((1.d0 + bhspin**2*rho2*hh)*psi4*rho2**2)/
     $    (rho2*bphi**2 +
     $    (alpha0**2 - bphi*bphiup)*(1.d0 + bhspin**2*rho2*
     $     hh)*psi4*rho2**2))

         gutx = -((y1*bphi*rho2)/(rho2*bphi**2 +
     $    (alpha0**2 - bphi*bphiup)*(1.d0 +
     $    bhspin**2*rho2*hh)*psi4*rho2**2))

         guty = (x1*bphi*rho2)/(rho2*bphi**2 +
     $    (alpha0**2 - bphi*bphiup)*(1.d0 +
     $    bhspin**2*rho2*hh)*psi4*rho2**2)

         gutz = 0.d0

         guxx = -((-(x1**2*bphi**2) -
     $    (alpha0**2 - bphi*bphiup)*(1.d0 +
     $    bhspin**2*x1**2*hh)*psi4*rho2**2)/(psi4*(rho2*
     $    bphi**2 + (alpha0**2 -
     $    bphi*bphiup)*(1.d0 + bhspin**2*rho2*hh)*
     $    psi4*rho2**2)))

         guyy = -((-(y1**2*bphi**2) -
     $    (alpha0**2 - bphi*bphiup)*(1.d0 +
     $    bhspin**2*y1**2*hh)*psi4*rho2**2)/(psi4*
     $    (rho2*bphi**2 + (alpha0**2 -
     $    bphi*bphiup)*(1.d0 + bhspin**2*rho2*hh)*psi4*rho2**2)))

         guzz = 1.d0 / psi4

         guxy =-((x1*y1*(-bphi**2 + bhspin**2*(-alpha0**2 +
     $   bphi*bphiup)*hh*psi4*rho2**2))/(psi4*
     $    (rho2*bphi**2 + (alpha0**2 -
     $   bphi*bphiup)*(1.d0 + bhspin**2*rho2*hh)*psi4*rho2**2)))

         guyz = 0.d0

         guzx = 0.d0




C-----------------------------------------------------------------------
C  Debug check: is gu^{mu nu} really the inverse of gd_{mu nu}?
C-----------------------------------------------------------------------

C        Fill covariant metric gdn_{mu nu}
         gdn(0,0) = gdtt
         gdn(0,1) = gdtx
         gdn(0,2) = gdty
         gdn(0,3) = gdtz
         gdn(1,0) = gdtx
         gdn(1,1) = gdxx
         gdn(1,2) = gdxy
         gdn(1,3) = gdzx
         gdn(2,0) = gdty
         gdn(2,1) = gdxy
         gdn(2,2) = gdyy
         gdn(2,3) = gdyz
         gdn(3,0) = gdtz
         gdn(3,1) = gdzx
         gdn(3,2) = gdyz
         gdn(3,3) = gdzz

C        Fill contravariant metric gu^{mu nu}
         gun(0,0) = gutt
         gun(0,1) = gutx
         gun(0,2) = guty
         gun(0,3) = gutz
         gun(1,0) = gutx
         gun(1,1) = guxx
         gun(1,2) = guxy
         gun(1,3) = guzx
         gun(2,0) = guty
         gun(2,1) = guxy
         gun(2,2) = guyy
         gun(2,3) = guyz
         gun(3,0) = gutz
         gun(3,1) = guzx
         gun(3,2) = guyz
         gun(3,3) = guzz

C        Compute product prod = gu * gd
         do i = 0, 3
            do j = 0, 3
               prod(i,j) = 0.d0
               do k = 0, 3
                  prod(i,j) = prod(i,j) + gun(i,k) * gdn(k,j)
               enddo
            enddo
         enddo

C        Check max deviation from identity
         tol    = 1.d-8
         maxerr = 0.d0
         do i = 0, 3
            do j = 0, 3
               if (i .eq. j) then
                  err = prod(i,j) - 1.d0
               else
                  err = prod(i,j)
               endif
               if (abs(err) .gt. maxerr) maxerr = abs(err)
            enddo
         enddo

         

         if (maxerr .gt. tol) then
            write(msg,'(A,3ES22.14,A,ES22.14)') 
     $   'Exact__Kerr_Quasi_isotropic: inverse check failed; x,y,z=',
     $        x, y, z, '  maxerr = ', maxerr
            call CCTK_WARN(0, msg)

            call CCTK_ERROR('Exact__Kerr_Quasi_isotropic: inverse metric check failed')
         endif
C-----------------------------------------------------------------------

         return
      end
