#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define Pi (4 * atan(1.d0))


      subroutine  gauge_wave(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx)

      implicit none
      integer CCTK_Equals

      DECLARE_CCTK_PARAMETERS

      CCTK_REAL x, y, z, t
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx

      logical firstcall

      CCTK_REAL a, o, H, d, fs
      CCTK_REAL zero,half,one

      data firstcall /.true./
      save firstcall, a, o, half, d

C     Numbers.

      zero = 0.0d0
      half = 0.5d0
      one  = 1.0d0

C     Get parameters of the exact solution.

      if (firstcall) then

         a = GW_a
         o = GW_omega
         d = GW_del
         fs = GW_phase_shift

         firstcall = .false.

      end if

C     How should the wave look like.

      if (CCTK_Equals(GW_H,"sin").eq.1) then

        d = GW_del * half / Pi

        if (GW_diagonal.ne.0) then

          H = one - a * sin((x-y)/sqrt(2.d0)/d - o*t/d)
 
        else

          H = one - a * sin((x-o*t)/d - fs)

        end if

      elseif (CCTK_Equals(GW_H,"exp").eq.1) then

        H = exp(a*(x/d-o*t))

      elseif (CCTK_Equals(GW_H,"gau").eq.1) then

        H = one - a*dexp(-(x-t)**2/d**2)

      end if

C     write metric.

      if (GW_diagonal .eq. 1) then

        gdxx = half * H + half
        gdxy = -half* H + half
        gdyy = half * H + half

        guxx = half / H + half
        guxy =-half / H + half
        guyy = half / H + half

      else

        gdxx = H
        gdxy = zero
        gdyy = one

        guxx = one/H
        guxy = zero
        guyy = one

      end if
      
      gdtt = - H
      gdtx = zero
      gdty = zero
      gdtz = zero

      gdzx = zero
      gdyz = zero
      gdzz = one

C     and upper metric.

      gutt = - one/H
      gutx = zero
      guty = zero
      gutz = zero

      guyz = zero
      guzx = zero
      guzz = one

      return
      end
