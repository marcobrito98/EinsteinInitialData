C     Kerr-Schild form of boosted rotating black hole.
C     Program g_ab = eta_ab + H l_a l_b, g^ab = eta^ab - H l^a l^b.
C     Here eta_ab is Minkowski in Cartesian coordinates, H is a scalar,
C     and l is a null vector.
C
C Author: unknown
C Copyright/License: unknown
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

      subroutine Exact__Kerr_KerrSchild(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

c input arguments
      CCTK_REAL x, y, z, t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_REAL psi
      LOGICAL   Tmunu_flag

c local static variables
      logical firstcall
      CCTK_REAL boostv, eps, m, a
      data firstcall /.true./
      save firstcall, boostv, eps, m, a

c local variables
      CCTK_REAL gamma, t0, z0, x0, y0, rho02, r02, r0, costheta0, 
     $     lt0, lx0, ly0, lz0, hh, lt, lx, ly, lz

      integer R02_TOO_SMALL_WARN_LEVEL
      parameter (R02_TOO_SMALL_WARN_LEVEL = 2)

      character*100 warn_buffer

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution,
C     and convert from parameter file spin parameter J/m^2
C     to the J/m definition used in the code here.

      if (firstcall) then
         boostv = Kerr_KerrSchild__boost_v
         eps    = Kerr_KerrSchild__epsilon
         m      = Kerr_KerrSchild__mass
         a      = m*Kerr_KerrSchild__spin
         firstcall = .false.
      end if

C     Boost factor.

      gamma = 1.d0 / sqrt(1.d0 - boostv**2)

C     Lorentz transform t,x,y,z -> t0,x0,y0,z0. 
C     t0 is never used, but is here for illustration, and we introduce
C     x0 and y0 also only for clarity.
C     Note that z0 = 0 means z = vt for the BH.

      t0 = gamma * (t - boostv * z)
      z0 = gamma * (z - boostv * t)
      x0 = x
      y0 = y

C     Coordinate distance to center of black hole. Note it moves!

      rho02 = x0**2 + y0**2 + z0**2

C     Spherical auxiliary coordinate r and angle theta in BH rest frame.

      r02 = 0.5d0 * (rho02 - a**2) 
     $     + sqrt(0.25d0 * (rho02 - a**2)**2 + a**2 * z0**2)
      if (r02 .lt. eps) then
         call CCTK_WARN(R02_TOO_SMALL_WARN_LEVEL,'Kerr/Kerr-Schild: r02 is too small at relative position')

         write(warn_buffer, '(a,f8.3, a,f8.3, a,f8.3)')
     $        '                  x0=',x0, ' y0=',y0, ' z0=',z0
         call CCTK_WARN(R02_TOO_SMALL_WARN_LEVEL, warn_buffer)

         write(warn_buffer, '(a,f9.3, a,g12.3)')
     $        '                  rho02=',rho02, ' ==> r02=',r02
         call CCTK_WARN(R02_TOO_SMALL_WARN_LEVEL, warn_buffer)
      end if
      r02 = max(r02,eps)
      r0 = sqrt(r02)
      costheta0 = z0 / r0
      
C     Coefficient H. Note this transforms as a scalar, so does not carry
C     the suffix 0.
      hh = m * r0 / (r0**2 + a**2 * costheta0**2)
      
C     Components of l_a in rest frame. Note indices down.
      lt0 = 1.d0
      lx0 = (r0 * x0 + a * y0) / (r0**2 + a**2)
      ly0 = (r0 * y0 - a * x0) / (r0**2 + a**2)
      lz0 = z0 / r0

C     Now boost it to coordinates x, y, z, t.
C     This is the reverse Lorentz transformation, but applied
C     to a one-form, so the sign of boostv is the same as the forward
C     Lorentz transformation applied to the coordinates.

      lt = gamma * (lt0 - boostv * lz0) 
      lz = gamma * (lz0 - boostv * lt0) 
      lx = lx0
      ly = ly0
      
C     Down metric. g_ab = flat_ab + H l_a l_b

      gdtt = - 1.d0 + 2.d0 * hh * lt * lt
      gdtx = 2.d0 * hh * lt * lx
      gdty = 2.d0 * hh * lt * ly
      gdtz = 2.d0 * hh * lt * lz
      gdxx = 1.d0 + 2.d0 * hh * lx * lx
      gdyy = 1.d0 + 2.d0 * hh * ly * ly
      gdzz = 1.d0 + 2.d0 * hh * lz * lz
      gdxy = 2.d0 * hh * lx * ly
      gdyz = 2.d0 * hh * ly * lz
      gdzx = 2.d0 * hh * lz * lx

C     Up metric. g^ab = flat^ab - H l^a l^b.
C     Notice that g^ab = g_ab and l^i = l_i and l^0 = - l_0 in flat spacetime.
      gutt = - 1.d0 - 2.d0 * hh * lt * lt 
      gutx = 2.d0 * hh * lt * lx
      guty = 2.d0 * hh * lt * ly
      gutz = 2.d0 * hh * lt * lz
      guxx = 1.d0 - 2.d0 * hh * lx * lx
      guyy = 1.d0 - 2.d0 * hh * ly * ly
      guzz = 1.d0 - 2.d0 * hh * lz * lz
      guxy = - 2.d0 * hh * lx * ly
      guyz = - 2.d0 * hh * ly * lz
      guzx = - 2.d0 * hh * lz * lx

      return
      end
