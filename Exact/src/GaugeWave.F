#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

C    #define Pi 4 * atan(1.d0)
C #define Pi 3.1415
      
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

      CCTK_REAL a, o, H, half, d, fs

      data firstcall /.true./
      save firstcall, a, o, half, d, fs

      

C     Get parameters of the exact solution.

      if (firstcall) then

         a = GW_a
         o = GW_omega
         half = 1/2d0
         d = GW_del / 8d0 / atan(1.d0)
         fs = GW_phase_shift

         firstcall = .false.

      end if

C     How should the wave look like
      if (CCTK_Equals(GW_H,"sin").eq.1) then
        
        if (GW_diagonal) then

          H = 1 - a * sin((x-y)/sqrt(2.d0)/d - o*t/d)
 
        else

          H = 1 - a * sin((x-o*t)/d - fs )

        end if

      elseif (CCTK_Equals(GW_H,"exp").eq.1) then

        H = exp(a*(x/d-o*t))

      end if

C     write metric      

      if (GW_diagonal) then
         gdxx = half + half *H 
         gdyy = half + half *H 
         gdxy = half - half *H  
         gdzz = 1
         gdzx = 0
         gdyz = 0

         gdtt = - H
         gdtx = 0
         gdty = 0
         gdtz = 0

C     and upper metric

         guxx = half + half / H 
         guyy = half + half / H
         guxy = half - half / H
         guzz = 1
         guyz = 0
         guzx = 0
      
         gutt = - 1/H
         gutx = 0
         guty = 0
         gutz = 0

      else

         gdxx = H 
         gdyy = 1
         gdzz = 1
         gdxy = 0 
         gdzx = 0
         gdyz = 0

         gdtt = - H
         gdtx = 0
         gdty = 0
         gdtz = 0

C     upper metric

         guxx = 1/H
         guyy = 1
         guzz = 1
         guxy = 0
         guyz = 0
         guzx = 0
              
         gutt = - 1/H
         gutx = 0
         guty = 0
         gutz = 0
      end if

      return
      end
