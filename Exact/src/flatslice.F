#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine flatslice(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx)

c     The metric given here corresponds to that of flat spacetime
c     but with a non-trivial slicing.  Basically, I take
c     the flat metric in spherical coordinates, and define
c     a new time coordinate such that:
c                                   
c     t  =  t    +  f(r)    =>   dt  =  dt    + f'(r)  dr
c            new                          new
c
c     Notice that only f'(r) ever appear in the new metric.  For
c     f'(r) I take a gaussian with amplitude a centered at r=r0.
c     For  -1 < a < 1, the transformation above gives us a
c     spacelike slicing.

      implicit none

      DECLARE_CCTK_PARAMETERS

c     Input.

      CCTK_REAL x,y,z,t

c     Output.

      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $          gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $          gutt, gutx, guty, gutz, 
     $          guxx, guyy, guzz, guxy, guyz, guzx

c     Internal.

      CCTK_REAL a,s,p,fp
      CCTK_REAL r,r2,x2,y2,z2
      CCTK_REAL det,aux
      CCTK_REAL zero,one,two

c     Define numbers.

      zero = 0.0d0
      one  = 1.0d0
      two  = 2.0d0

c     Read parameters

      a = flatslice_a
      s = flatslice_s
      p = flatslice_r

c     Find transformation function.  Notice that I use
c     a gaussian times r^2/(1+r^2).  The extra factor
c     is there to make sure that the transformation is
c     regular at the origin.

      x2 = x**2
      y2 = y**2
      z2 = z**2

      r2 = max(x2 + y2 + z2,1.0d-20)
      r  = dsqrt(r2)

c      fp = a*r2/(one+r2)*dexp(-(r-p)**2/s**2)
      fp = a*r2/(p+r2)*exp(-r/s)

c     Give metric components.

      gdtt = - one

      gdtx = - fp*x/r
      gdty = - fp*y/r
      gdtz = - fp*z/r

      aux = fp**2/r2

      gdxx = one - x2*aux
      gdyy = one - y2*aux
      gdzz = one - z2*aux

      gdxy = - x*y*aux
      gdyz = - y*z*aux
      gdzx = - z*x*aux

c     Find inverse metric.

      gutt = - one + fp**2

      gutx = - fp*x/r
      guty = - fp*y/r
      gutz = - fp*z/r

      det = gdxx*gdyy*gdzz + two*gdxy*gdzx*gdyz
     .    - gdxx*gdyz**2 - gdyy*gdzx**2 - gdzz*gdxy**2

      guxx = (gdyy*gdzz - gdyz**2)/det
      guyy = (gdxx*gdzz - gdzx**2)/det
      guzz = (gdxx*gdyy - gdxy**2)/det

      guxy = (gdzx*gdyz - gdzz*gdxy)/det
      guyz = (gdxy*gdzx - gdxx*gdyz)/det
      guzx = (gdxy*gdyz - gdyy*gdzx)/det

      return
      end
