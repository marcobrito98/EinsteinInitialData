C DeSitter metric spacetime with cosmological constant
C Author : D. Vulcanov (Timisoara, Romania)
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__de_Sitter_Lambda(
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
    
      CCTK_REAL arad,am

      data firstcall /.true./
      save firstcall, arad
      if (firstcall) then
           arad = de_Sitter_Lambda__scale
          firstcall = .false.
      end if

      am=exp(2.0D0*dsqrt(arad/3.0d0)*t)   

      gdtt = -1.d0 
      gdtx = 0.d0 
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = am
      gdyy = am
      gdzz = am
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -1.d0
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = 1.d0/am
      guyy = 1.d0/am
      guzz = 1.d0/am
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0


      return
      end
