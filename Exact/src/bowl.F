#include "cctk.h"
#include "cctk_parameters.h"

      subroutine bowl(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx)

c     The metric given here is not a solution of
c     Einstein's equations!  It is nevertheless
c     useful for tests since it has a particularly
c     nice geometry.  It is a static, spherically
c     symmetric metric with no shift.
c
c     In spherical coordinates, the metric has the
c     form:
c
c       2      2        2        2
c     ds  =  dr  +  R(r)  d Omega
c
c     Clearly, r measures radial proper distance, and R(r)
c     is the areal (Schwarzschild) radius.
c
c     I choose a form of R(r) such that:
c
c     R  -->   r     r<<1, r>>1
c
c     So close in, and far away we have a flat metric.
c     In the middle region, I take R to be smaller than
c     r, but still larger than zero.  This deficit in
c     areal radius produces the geometry of a `bag of gold'.
c
c     The size of the deviation from a flat geometry
c     is controled by the parameter "bowl_a".
c     If this parameter is 0, we are in flat space.
c     The width of the curved region is controled by
c     the paramter "bowl_s", and the place where the
c     curvature becomes significant (the center of the
c     deformation) is controled by "bowl_c".
c
c     The specific form of the function R(r) is:
c
c     R(r) = ( r - a f(r) )
c
c     and the form of thr function f(r) depends on the
c     parameter bowl_type:
c                                            2   2
c     bowl_type = "gauss"   f(r) = exp(-(r-c) / s )
c
c     bowl_type = "fermi"   f(r) = 1 / ( 1 + exp(-s(r-c)) )
c
c     There are three extra paramters (bowl_dx,bowl_dy,bowl_dz)
c     that scale the (x,y,z) axis respectively.  Their default
c     value is 1.  These parameters are useful to hide the
c     spherical symmetry of the metric.

      implicit none

      DECLARE_CCTK_PARAMETERS

c     Input.
      CCTK_REAL x, y, z, t

c     Output.
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $          gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $          gutt, gutx, guty, gutz, 
     $          guxx, guyy, guzz, guxy, guyz, guzx

      logical firstcall,evolve

      integer type
      integer CCTK_Equals

      CCTK_REAL a,c,s
      CCTK_REAL dx,dy,dz
      CCTK_REAL r,r2,rr2
      CCTK_REAL xr,yr,zr,xr2,yr2,zr2
      CCTK_REAL fac,det
      CCTK_REAL tfac,st,t0
      CCTK_REAL zero,one,two

      data firstcall /.true./
      save firstcall,evolve,type,a,c,s,dx,dy,dz,t0,st

c     Find {zero,one,two}.

      zero = 0.0d0
      one  = 1.0d0
      two  = 2.0d0

c     Get parameters of the metric.

      if (firstcall) then

         a = bowl_a
         c = bowl_c
         s = bowl_s

         dx = bowl_dx
         dy = bowl_dy
         dz = bowl_dz

         print *, a,c,s,dx,dy,dz

         if (dx.le.zero) then
            call CCTK_WARN(0,"The parameter bowl_dx must be greater than 0")
         end if

         if (dy.le.zero) then
            call CCTK_WARN(0,"The parameter bowl_dy must be greater than 0")
         end if

         if (dz.le.zero) then
            call CCTK_WARN(0,"The parameter bowl_dz must be greater than 0")
         end if

         if (CCTK_Equals(bowl_type,"gauss").ne.0) then
            type = 1
         else if (CCTK_Equals(bowl_type,"fermi").ne.0) then
            type = 2
         else
            call CCTK_WARN(0,"Unknown type of bowl metric")
         end if

         if (bowl_evolve.eq.1) then
            evolve = .true.
            t0 = bowl_t0
            st = bowl_st
         else
            evolve = .false.
         end if

         firstcall = .false.

      end if

c     Multiply the bowl strength `a' with a time evolution factor.
c     The time evolution factor is taken to be a Fermi step centered
c     in `t0' and with a width `st'.  The size of this step is always
c     1 so that far in the past we will always have flat space, and
c     far in the future we will have the static bowl.

      if (evolve) then
         tfac = one/(one + dexp(-st*(t-t0)))
      else
         tfac = one
      end if

      a = a*tfac

c     Find {r2,r}.

      r2 = (x/dx)**2 + (y/dy)**2 + (z/dz)**2
      r  = dsqrt(r2)

c     Find the form function rr2
c
c                      2    2                 2
c     rr2  =  (r - a f)  / r  =  (1 - a f / r)

      if (type.eq.1) then

c        Gaussian bowl:
c                                2  2       2
c        rr2 = [ 1 - a exp(-(r-c) /s ) / r ]
c
c        Notice that this really does not go to 1 at the
c        origin.  To fix this, I multiply the gaussian
c        with the factor:
c
c        fac  =  1  -  sech(4r)
c
c        This goes smoothly to 0 at the origin, and climbs
c        fast to a limiting value of 1 (at r=1 it is already
c        equal to 0.96).

         fac = one - two/(dexp(4.0d0*r) + dexp(-4.0d0*r))
         rr2 = (one - a*fac*dexp(-((r-c)/s)**2)/r)**2

      else if (type.eq.2) then

c        Fermi bowl:
c                                                  2
c        rr2 = [ 1 - 1 / ( 1 + exp(-s(r-c)) ) / r ]
c
c        Again, this doesn't really go to 1 at the origin, so
c        I use the same trick as above.

         fac = one - two/(dexp(4.0d0*r) + dexp(-4.0d0*r))
         rr2 = (one - a*fac/(one + dexp(-s*(r-c)))/r)**2

      end if

c     Give metric components.

      gdtt = - one
      gdtx = zero
      gdty = zero
      gdtz = zero

      if (r.ne.0) then

         xr = (x/dx)/r
         yr = (y/dy)/r
         zr = (z/dz)/r

         xr2 = xr**2
         yr2 = yr**2
         zr2 = zr**2

         gdxx = (xr2 + rr2*(yr2 + zr2))/dx**2
         gdyy = (yr2 + rr2*(xr2 + zr2))/dy**2
         gdzz = (zr2 + rr2*(xr2 + yr2))/dz**2

         gdxy = xr*yr*(one - rr2)/(dx*dy)
         gdyz = yr*zr*(one - rr2)/(dy*dz)
         gdzx = xr*zr*(one - rr2)/(dx*dz)

      else

         gdxx = one
         gdyy = one
         gdzz = one

         gdxy = zero
         gdyz = zero
         gdzx = zero

      end if

c     Find inverse metric.

      gutt = - one
      gutx = zero
      guty = zero
      gutz = zero

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
