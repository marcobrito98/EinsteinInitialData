C Einstein-DeSitter metric spacetime !!!!
C It emulates the Robertson-Walker universe
C near t=0, with zero pressure, and k=0
C See :J.N. Islam, An Introduction to 
C Mathematical Cosmology, Cambridge, 1992 and 
C S. Hawking, G.F.R. Ellis, The Large Scale
C Structure of space-time, Cambridge, 1973
C Author : D. Vulcanov (Timisoara, Romania)

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine  desitter(
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
    
      CCTK_REAL arad,r2,x2,y2,z2,am

      data firstcall /.true./
      save firstcall
      if (firstcall) then
          arad = desitt_a
          PRINT *, firstcall, arad, desitt_a
          firstcall = .false.
      end if
 
      x2=x*x
      y2=y*y
      z2=z*z
      r2 = x2+y2+z2
      am=arad*t**(4/3)     

      gdtt = -1.d0 
      gdtx = 0.d0 
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = 1.d0+(-1.d0 +am)*x2/r2 
      gdyy = 1.d0+(-1.d0 +am)*y2/r2
      gdzz = 1.d0+(-1.d0 +am)*z2/r2
      gdxy = (-1.d0   +am)*x*y/r2
      gdyz = (-1.d0   +am)*y*z/r2
      gdzx = (-1.d0   +am)*x*z/r2

      gutt = -1.d0
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = 1.d0+(-1.d0 +1.d0/am)*x2/r2 
      guyy = 1.d0+(-1.d0 +1.d0/am)*y2/r2
      guzz = 1.d0+(-1.d0 +1.d0/am)*z2/r2
      guxy = (-1.d0  +1.d0/am)*x*y/r2
      guyz = (-1.d0  +1.d0/am)*y*z/r2
      guzx = (-1.d0  +1.d0/am)*x*z/r2
 


      return
      end
