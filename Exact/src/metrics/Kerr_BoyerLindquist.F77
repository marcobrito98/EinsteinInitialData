C Kerr metric in cartesian Boyer-Lindquist coordinates,
C as per MTW box 33.2.
C Author : D. Vulcanov (Timisoara, Romania)
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Kerr_BoyerLindquist(
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
    
      CCTK_REAL arad, marad

      data firstcall /.true./
      save firstcall, arad, marad
      if (firstcall) then
          arad = Kerr_BoyerLindquist__spin
         marad = Kerr_BoyerLindquist__mass
          firstcall = .false.
      end if


     
      gdtt = -(y**2*arad**2+x**2-2*marad*x)/(x**2+y**2*arad**2)
      gdtx = 2*(arad*marad*x*(y**2-1))/(x**2+y**2*arad**2)
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = -(x**4+x**2*arad**2+2*arad**2*marad*x+arad**2*y**2*x**2 - 2*arad**2*y**2*marad*x+arad**4*y**2)*(y**2-1)/(x**2+arad**2*y**2)   
      gdyy = (x**2+y**2*arad**2)/(x**2-2*marad*x+arad**2) 
      gdzz =  -(x**2+y**2*arad**2)/(y**2-1)
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -(-x**4-x**2*arad**2-2*arad**2*marad*x-arad**2*y**2*x**2 +2*arad**2*y**2*marad*x-arad**4*y**2)/((x**2+arad**2*y**2)*(-x**2+2*marad*x-arad**2))
      gutx = 2*(arad*marad*x)/((x**2+arad**2*y**2)*(-x**2+2*marad*x-arad**2))
      guty = 0.d0
      gutz = 0.d0
      guxx = -(-arad**2*y**2-x**2+2*marad*x)/((x**2+arad**2*y**2)*(y**2-1)*(-x**2+2*marad*x-arad**2))
      guyy = -(-x**2+2*marad*x-arad**2)/(x**2+arad**2*y**2)
      guzz = -(y**2-1)/(x**2+arad**2*y**2)
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0

      return
      end
