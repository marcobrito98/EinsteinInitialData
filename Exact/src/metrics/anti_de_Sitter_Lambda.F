C Anti DeSitter metric spacetime with csomological constant
C Author : D. Vulcanov (Timisoara, Romania)
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__anti_de_Sitter_Lambda(
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
          arad = anti_de_Sitter_Lambda__scale
          firstcall = .false.
      end if
 


      am=exp(2.0D0*dsqrt(-arad/3.0D0)*x)   

       
      gdtt = -am 
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = 1.0D0
      gdyy = am
      gdzz = am
      gdxy = 0.0D0
      gdyz = 0.0D0
      gdzx = 0.0D0

      gutt = -1.0D0/am
      gutx = 0.0D0
      guty = 0.0D0
      gutz = 0.0D0
      guxx = 1.0D0
      guyy = 1.0D0/am
      guzz = 1.0D0/am
      guxy = 0.0D0
      guyz = 0.0D0
      guzx = 0.0D0


      return
      end
