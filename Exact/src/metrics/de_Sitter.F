C Einstein-DeSitter metric spacetime !!!!
C It emulates the Robertson-Walker universe
C near t=0, with zero pressure, and k=0
C See :J.N. Islam, An Introduction to 
C Mathematical Cosmology, Cambridge, 1992 and 
C S. Hawking, G.F.R. Ellis, The Large Scale
C Structure of space-time, Cambridge, 1973
C Author : D. Vulcanov (Timisoara, Romania)
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__de_Sitter(
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
          arad = de_Sitter__scale
          firstcall = .false.
      end if
 
      am=arad*t**(4.0D0/3.0D0)     

      gdtt = -1.0D0 
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
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
