C Axisymmetric Kasner solution !
C Author : D. Vulcanov
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Kasner_axisymmetric(
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

      gdtt = -t**(-0.5d0) 
      gdtx = 0.d0 
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = t**(-0.5d0)
      gdyy = t
      gdzz = t
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -t**(0.5d0)
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = t**(0.5d0)
      guyy = 1.d0/t
      guzz = 1.d0/t
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0
 


      return
      end
