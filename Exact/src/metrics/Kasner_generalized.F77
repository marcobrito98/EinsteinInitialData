C Kasner generalized metric spacetime !!!!
C Author : D. Vulcanov
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Kasner_generalized(
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
    
      CCTK_REAL pp1,pp2,a1,a2,a3

      data firstcall /.true./
      save firstcall, pp1, pp2
      if (firstcall) then
          pp1= Kasner_generalized__p1
          pp2= Kasner_generalized__p2
          firstcall = .false.
      end if
 
      a1= t**(2.d0*pp1)
      a2= t**(2.d0*pp2)
      a3= t**(2.d0-2.d0*pp1-2.d0*pp2)    


      gdtt = -1.d0 
      gdtx = 0.d0 
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = a1 
      gdyy = a2
      gdzz = a3
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -1.d0
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = 1.d0/a1
      guyy = 1.d0/a2
      guzz = 1.d0/a3
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0
 

      return
      end
