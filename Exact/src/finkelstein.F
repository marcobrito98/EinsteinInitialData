#include "cctk.h"
#include "cctk_parameters.h"

      subroutine  finkelstein(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx)

      implicit none

      DECLARE_CCTK_PARAMETERS

      CCTK_REAL x, y, z, t
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx

      logical firstcall

      CCTK_REAL eps, m

      data firstcall /.true./
      save firstcall, eps, m

      CCTK_REAL r

C     Get parameters of the exact solution.

      if (firstcall) then

         eps = kerrschild_eps
         m   = kerrschild_m

         firstcall = .false.

      end if

      r = max(sqrt(x**2 + y**2 + z**2), eps)

      gdtt = - (1.d0 - 2.d0 * m / r)
      gdtx = 2.d0 * m * x / r**2
      gdty = 2.d0 * m * y / r**2
      gdtz = 2.d0 * m * z / r**2
      gdxx = 1.d0 + 2.d0 * m * x**2 / r**3
      gdyy = 1.d0 + 2.d0 * m * y**2 / r**3
      gdzz = 1.d0 + 2.d0 * m * z**2 / r**3
      gdxy = 2.d0 * m * x * y / r**3
      gdyz = 2.d0 * m * y * z / r**3
      gdzx = 2.d0 * m * z * x / r**3

      gutt = - (1.d0 + 2.d0 * m / r)
      gutx = 2.d0 * m * x / r**2
      guty = 2.d0 * m * y / r**2
      gutz = 2.d0 * m * z / r**2
      guxx = 1.d0 - 2.d0 * m * x**2 / r**3
      guyy = 1.d0 - 2.d0 * m * y**2 / r**3
      guzz = 1.d0 - 2.d0 * m * z**2 / r**3
      guxy = - 2.d0 * m * x * y / r**3
      guyz = - 2.d0 * m * y * z / r**3
      guzx = - 2.d0 * m * z * x / r**3

      return
      end
