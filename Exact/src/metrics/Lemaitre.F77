C Lemaitre type universe - FRW with k =0, p=k rho and 
C cosmological constant
C $Header$


#include "cctk.h"
#include "cctk_Parameters.h"


      subroutine Exact__Lemaitre(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi)

      implicit none

      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS      


      CCTK_REAL x, y, z, t
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx,
     $       psi
      

      logical firstcall
    
      CCTK_REAL la, e0, ka, unu, ra, ra2, r0, Pii

      data firstcall /.true./
      save firstcall, ka, la, e0, r0
      if (firstcall) then
          ka = Lemaitre__kappa
          e0 = Lemaitre__epsilon0
          la = Lemaitre__Lambda 
          r0 = Lemaitre__R0
          firstcall = .false.
      end if

       Pii = acos(-1.0D0)

       unu = dsqrt(3.0D0*la)*t*(ka+1.0D0)/(2.0D0)

       ra  = r0*(cosh(unu)+dsqrt(1.0D0+8.0D0*Pii*e0/la)*sinh(unu))**
     &       (2.0D0/(3.0D0*ka+3.0D0)) 

       ra2 = ra*ra  
 
      gdtt = -1.0D0 
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = ra2
      gdyy = ra2
      gdzz = ra2
      gdxy = 0.0D0
      gdyz = 0.0D0
      gdzx = 0.0D0

      gutt = -1.0D0
      gutx = 0.0D0
      guty = 0.0D0
      gutz = 0.0D0
      guxx = 1.0D0/ra2
      guyy = 1.0D0/ra2
      guzz = 1.0D0/ra2
      guxy = 0.0D0
      guyz = 0.0D0
      guzx = 0.0D0
 


      return
      end
