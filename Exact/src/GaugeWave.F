#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

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

      CCTK_REAL a, o, H, half, d

      data firstcall /.true./
      save firstcall, a, o, half, d

      

C     Get parameters of the exact solution.

      if (firstcall) then

         a = GW_a
         o = GW_omega
         half = 1/2d0
         d = GW_del

         firstcall = .false.

      end if

C     How should the wave look like
      if (CCTK_Equals(GW_H,"sin").eq.1) then
        
        H = 1-a*sin((x-t)/d)

      elseif (CCTK_Equals(GW_H,"exp").eq.1) then

        H = exp(a*(x/d-o*t))

      end if

C     write metric      
      if (GW_diagonal) then

        gdxx = half * H + half
        gdxy = -half* H + half
        gdyy = half * H + half

        guxx = half / H + half
        guxy =-half / H + half
        guyy = half / H + half

      else

        gdxx = H
        gdxy = 0
        gdyy = 1

        guxx = 1 / H
        guxy = 0 
        guyy = 1

      end if

      
      gdtt = - H
      gdtx = 0
      gdty = 0
      gdtz = 0

      gdzx = 0
      gdyz = 0
      gdzz = 1
      
      
C     and upper metric      
      gutt = - 1/H
      gutx = 0
      guty = 0
      gutz = 0

      guyz = 0
      guzx = 0
      guzz = 1

      return
      end
