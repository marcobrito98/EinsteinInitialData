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

         r2 = max(x1**2 + y1**2 + z1**2, eps)
         r = sqrt(r2)
         rho2 = max(x1**2 + y1**2, eps)
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

         gutt = -(((1.d0 + bhspin**2*(x1**2 + y1**2)*hh)*psi4*rho2**2) /
     $    ((x1**2 + y1**2)*bphi**2 +
     $    (alpha0**2 - bphi*bphiup)*(1 + bhspin**2*(x1**2 +
     $    y1**2)*hh)*psi4*rho2**2))
         gutx = -((y1*bphi*rho2)/((x1**2 + y1**2)*bphi**2 + (alpha0**2 -
     $    bphi*bphiup)*(1 + bhspin**2*(x1**2 + y1**2)*hh)*psi4*rho2**2))
         guty = (x1*bphi*rho2)/((x1**2 + y1**2)*bphi**2 + (alpha0**2 - 
     $    bphi*bphiup)*(1.d0 +
     $    bhspin**2*(x1**2 + y1**2)*hh)*psi4*rho2**2)
         gutz = 0.d0
         guxx = -((-(x1**2*bphi**2) - (alpha0**2 - bphi*bphiup)*(1.d0 +
     $    bhspin**2*x1**2*hh)*psi4*rho2**2) /
     $    (psi4*((x1**2 + y1**2)*bphi**2 +
     $    (alpha0**2 - bphi*bphiup)*(1.d0 + bhspin**2*(x1**2 +
     $    y1**2)*hh)*psi4*rho2**2)))
         guyy = -((-(y1**2*bphi**2) - (alpha0**2 - bphi*bphiup)*(1.d0 +
     $    bhspin**2*y1**2*hh)*psi4*rho2**2) /
     $    (psi4*((x1**2 + y1**2)*bphi**2 +
     $    (alpha0**2 - bphi*bphiup)*(1.d0 + bhspin**2*(x1**2 +
     $    y1**2)*hh)*psi4*rho2**2)))
         guzz = 1.d0 / psi4
         guxy =-((x1*y1*(-bphi**2 + bhspin**2*(-alpha0**2 +
     $    bphi*bphiup)*hh*psi4*rho2**2))/(psi4*((x1**2 + y1**2)*bphi**2 +
     $    (alpha0**2 - bphi*bphiup)*(1 + bhspin**2*(x1**2 +
     $    y1**2)*hh)*psi4*rho2**2)))
         guyz = 0.d0
         guzx = 0.d0

         return
      end
