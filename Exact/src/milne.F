C De Milne spacetime metric 
C Suggested by Matteo Rossi and E. Onofri (Univ. di Parma, Italy)
C They inted to use thsi metric for simulating an Pre-Big-Bang
C Cosmology, as proposed by Veneziano some year ago 
C Author : D. Vulcanov (Timsoara, Romania)
#include "cctk.h"

      subroutine milne(
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
     $       guxx, guyy, guzz, guxy, guyz, guzx,
     $       coef


      coef= t**2/(1+ x**2 +y**2+ z**2)

      gdtt = -1.d0 
      gdtx = 0.d0  
      gdty = 0.d0 
      gdtz = 0.d0 
      gdxx =  coef*(1+y**2+z**2) 
      gdyy =  coef*(1+x**2+z**2) 
      gdzz =  coef*(1+x**2+y**2)
      gdxy = -coef*x*y 
      gdyz = -coef*y*z 
      gdzx = -coef*x*z

      gutt = -1.d0 
      gutx = 0.d0 
      guty = 0.d0 
      gutz = 0.d0 
      guxx = (1+x**2)/t**2
      guyy = (1+y**2)/t**2
      guzz = (1+z**2)/t**2
      guxy = x*y/t**2
      guyz = y*z/t**2 
      guzx = x*z/t**2


      return
      end
