C Bertotti spacetime with cosmological constant
C Author : D. Vulcanov (Timisoara, Romania)
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Bertotti(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi)

      implicit none

      DECLARE_CCTK_PARAMETERS

      CCTK_REAL x, y, z, t
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx,
     $       psi

      logical firstcall
    
      CCTK_REAL baza, unu, doi

      data firstcall /.true./
      save firstcall, baza
      if (firstcall) then
          baza = Bertotti__Lambda
          firstcall = .false.
      end if
 

      unu=exp(2.0D0*dsqrt(-baza)*x)
      doi=exp(2.0D0*dsqrt(-baza)*z)   

       
      gdtt = -unu
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = 1.0D0
      gdyy = doi
      gdzz = 1.0D0
      gdxy = 0.0D0
      gdyz = 0.0D0
      gdzx = 0.0D0

      gutt = -1.0D0/unu
      gutx = 0.0D0
      guty = 0.0D0
      gutz = 0.0D0
      guxx = 1.0D0
      guyy = 1.0D0/doi
      guzz = 1.0D0
      guxy = 0.0D0
      guyz = 0.0D0
      guzx = 0.0D0


      return
      end
