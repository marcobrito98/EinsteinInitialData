c $Header$

C     Author: Frank Loeffler (frank.loeffler@aei.mpg.de)
C     Licence: GPL 2 or later

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define Pi (4 * atan(1.d0))


      subroutine Exact__Minkowski_conf_wave(
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

C     write conformal factor
      if (Minkowski_conf_wave__direction .eq. 0) then
        psi = Minkowski_conf_wave__amplitude *
     .        sin(2.0d0*Pi/Minkowski_conf_wave__wavelength* x ) + 1.0d0
      else if (Minkowski_conf_wave__direction .eq. 1) then
        psi = Minkowski_conf_wave__amplitude *
     .        sin(2.0d0*Pi/Minkowski_conf_wave__wavelength* y ) + 1.0d0
      else if (Minkowski_conf_wave__direction .eq. 2) then
        psi = Minkowski_conf_wave__amplitude *
     .        sin(2.0d0*Pi/Minkowski_conf_wave__wavelength* z ) + 1.0d0
      end if

C     write metric.

      gdxx = psi**(-4.0d0)
      gdyy = gdxx
      gdzz = gdxx

      gdxy = 0.0d0
      gdyz = 0.0d0
      gdzx = 0.0d0

      gdtt = -1.0d0
      gdtx = 0.0d0
      gdty = 0.0d0
      gdtz = 0.0d0

C     and upper metric.

      guxx = psi**4.0d0
      guyy = guxx
      guzz = guxx

      guxy = 0.0d0
      guyz = 0.0d0
      guzx = 0.0d0

      gutt = -1.0d0
      gutx = 0.0d0
      guty = 0.0d0
      gutz = 0.0d0

      return
      end
