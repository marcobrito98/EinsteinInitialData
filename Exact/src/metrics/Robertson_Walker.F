C Robertson-Walker universe 
C See: J.N. Islam, An Introduction to 
C Mathematical Cosmology, Cambridge, 1992
C Author : D. Vulcanov (Timisoara, Romania)
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"


      subroutine Exact__Robertson_Walker(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,rama)

      implicit none

      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS      


      CCTK_REAL x, y, z, t
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx,rama
      

      logical firstcall
    
      CCTK_REAL r2,x2,y2,z2,am,ag, kapa

      data firstcall /.true./
      save firstcall, kapa
      if (firstcall) then
          kapa = Robertson_Walker__k
          firstcall = .false.
      end if

      x2=x*x
      y2=y*y
      z2=z*z
      r2 = x2+y2+z2

      ag = kapa/(1.0D0-kapa*r2)    
       
      am = rama*rama
       


      gdtt = -1.0D0 
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = am*(1.0D0 + ag*x2) 
      gdyy = am*(1.0D0 + ag*y2)
      gdzz = am*(1.0D0 + ag*z2)
      gdxy = am*ag*x*y
      gdyz = am*ag*y*z
      gdzx = am*ag*z*x

      gutt = -1.0D0
      gutx = 0.0D0
      guty = 0.0D0
      gutz = 0.0D0
      guxx = (1.0D0-kapa*x2)/am 
      guyy = (1.0D0-kapa*y2)/am
      guzz = (1.0D0-kapa*z2)/am
      guxy = -kapa*x*y/am
      guyz = -kapa*y*z/am
      guzx = -kapa*z*x/am
 


      return
      end
