C Kasner_like metric spacetime !!!!
C See: L. Pimentel, Int. Journ. of Theor. Physics,
C, vol. 32, No. 6, 1993, p. 979 (and the references
C cited here), and S. Gotlober, et. al. Early Evolution
C of the Universe and Formation Structure, Akad. Verlag, 1990
C Author : D. Vulcanov

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine  kasner(
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
    
      CCTK_REAL qq,a1,a3

      data firstcall /.true./
      save firstcall
      if (firstcall) then
          qq = kasner_q
          PRINT *, firstcall, qq, kasner_q
          firstcall = .false.
      end if
 
      a1= t**qq
      a3= t**(1-2*qq)    

      gdtt = -1.d0 
      gdtx = 0.d0 
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = a1**2 
      gdyy = a1**2
      gdzz = a3**2
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -1.d0
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = a1**(-2)
      guyy = a1**(-2)
      guzz = a3**(-2)
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0
 


      return
      end
