#include "cctk.h"
#include "cctk_parameters.h"

      subroutine starschwarz(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx)

c     The metric given here corresponds to a constant
c     density star, also known as a "Schwarzschild" star.
c     It is obviously not a vacuum solution, so if someone
c     wants to evolve it, we will need to add here the
c     matter variables.  For the moment, only the metric
c     is given.
c
c     The metric is given as a conformally flat metric.
c     Turns out that in the original areal radius, the
c     metric variables have a kink at the surface of the
c     star, but they are smooth in the conformal form.
c
c     Thanks to Philippos Papadopoulos for suggesting
c     the use of this metric.

      implicit none

      DECLARE_CCTK_PARAMETERS

c     Input.

      CCTK_REAL x, y, z, t

c     Output.

      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx

c     Internal.

      logical firstcall

      CCTK_REAL mass,radius
      CCTK_REAL r,c,psi4
      CCTK_REAL zero,one,two

      data firstcall /.true./
      save firstcall

c     Find {zero,one,two}.

      zero = 0.0D0
      one  = 1.0D0
      two  = 2.0D0

c     Get parameters of the metric.

      if (firstcall) then

         mass   = starschwarz_m
         radius = starschwarz_r

         firstcall = .false.

      end if

c     Find r.

      r = dsqrt(x**2 + y**2 + z**2)

c     Find conformal factor.

      if (r.le.radius) then

         c = mass/(two*radius)

         psi4 = (one + c)**6/(one + c*(r/radius)**2)**2

      else

         c = mass/(two*r)

         psi4 = (one + c)**4

      end if

c     Find metric components.

      gdtt = - psi4

      gdtx = zero
      gdty = zero
      gdtz = zero

      gdxx = psi4
      gdyy = psi4
      gdzz = psi4

      gdxy = zero
      gdyz = zero
      gdzx = zero

c     Find inverse metric.

      gutx = zero
      guty = zero
      gutz = zero

      guxx = one/psi4
      guyy = one/psi4
      guzz = one/psi4

      guxy = zero
      guyz = zero
      guzx = zero

      return
      end
