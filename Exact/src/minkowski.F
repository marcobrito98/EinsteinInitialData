#include "cctk.h"

      subroutine minkowski(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx)

      implicit none

c     Input.

      CCTK_REAL x, y, z, t

c     Output.

      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx

      gdtt = -1.d0 
      gdtx = 0.d0 
      gdty = 0.d0 
      gdtz = 0.d0 
      gdxx = 1.d0 
      gdyy = 1.d0 
      gdzz = 1.d0 
      gdxy = 0.d0 
      gdyz = 0.d0 
      gdzx = 0.d0

      gutt = -1.d0 
      gutx = 0.d0 
      guty = 0.d0 
      gutz = 0.d0 
      guxx = 1.d0 
      guyy = 1.d0 
      guzz = 1.d0 
      guxy = 0.d0 
      guyz = 0.d0 
      guzx = 0.d0

      return
      end
