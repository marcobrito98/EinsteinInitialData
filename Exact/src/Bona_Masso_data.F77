c This routine calculates Bona-Masso initial data, making use of the
c subroutine Exact__metric() to calculate the spacetime metric and its
c inverse.  Note that this use of the Bona-Masso variables is independent
c of how (or even if) we are evolving the Einstein equations -- here
c the Bona-Masso variables are "just" used as intermediate variables.

c $Header$

#include "cctk.h"

      subroutine  Exact__Bona_Masso_data(
     $     decoded_exact_model,
     $     x, y, z, t,
     $     gxx, gyy, gzz, gxy, gyz, gxz,
     $     hxx, hyy, hzz, hxy, hyz, hxz,
     $     dxgxx, dxgyy, dxgzz, dxgxy, dxgyz, dxgxz,
     $     dygxx, dygyy, dygzz, dygxy, dygyz, dygxz,
     $     dzgxx, dzgyy, dzgzz, dzgxy, dzgyz, dzgxz,
     $     alp, ax, ay, az, betax, betay, betaz,
     $     bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz)

      implicit none
      CCTK_INT decoded_exact_model
      CCTK_REAL x, y, z, t,
     $       gxx, gyy, gzz, gxy, gyz, gxz,
     $       hxx, hyy, hzz, hxy, hyz, hxz,
     $       dxgxx, dxgyy, dxgzz, dxgxy, dxgyz, dxgxz,
     $       dygxx, dygyy, dygzz, dygxy, dygyz, dygxz,
     $       dzgxx, dzgyy, dzgzz, dzgxy, dzgyz, dzgxz,
     $       alp, ax, ay, az, betax, betay, betaz,
     $       bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz

C     gxx is g_xx etc.
C     hxx is K_xx etc.
C     dxgyy is (/2) dg_yy / dx (sic!)
C     alp is N, betax is N^x etc.
C     bxy is (/2) dN^y / dx (sic and sic!)
C     ax is dN / dx / N (sic!)

      CCTK_REAL eps, 
     $       gdtt, gdtx, gdty, gdtz, 
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guxz,
     $       gdtt_p, gdtx_p, gdty_p, gdtz_p, 
     $       gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $       gutt_p, gutx_p, guty_p, gutz_p, 
     $       guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $       gdtt_m, gdtx_m, gdty_m, gdtz_m, 
     $       gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $       gutt_m, gutx_m, guty_m, gutz_m, 
     $       guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m

      parameter (eps=1.d-6)

C     Get the spacetime metric and its inverse at the base point.

      call Exact__metric(
     $     decoded_exact_model,
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gxx, gyy, gzz, gxy, gyz, gxz,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guxz)

C     Calculate lapse and shift from the upper metric.

      alp = 1.d0 / sqrt(-gutt)

      betax = - gutx / gutt
      betay = - guty / gutt
      betaz = - gutz / gutt

C     In order to calculate the derivatives of the lapse and shift from
C     the contravariant metric, use g^tt = -1/N^2 and g^i0 = N^i / N^2
C     Note that ax is dN/dx / N and that bxy is 1/2 dN^y / dx.

C     Calculate x-derivatives. 

      call Exact__metric(
     $     decoded_exact_model,
     $     x+eps, y, z, t,
     $     gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $     gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $     gutt_p, gutx_p, guty_p, gutz_p,
     $     guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p)
      call Exact__metric(
     $     decoded_exact_model,
     $     x-eps, y, z, t,
     $     gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $     gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $     gutt_m, gutx_m, guty_m, gutz_m,
     $     guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m)

      dxgxx = (gdxx_p  - gdxx_m)  / 4.d0 / eps
      dxgyy = (gdyy_p  - gdyy_m)  / 4.d0 / eps
      dxgzz = (gdzz_p  - gdzz_m)  / 4.d0 / eps
      dxgxy = (gdxy_p  - gdxy_m)  / 4.d0 / eps
      dxgyz = (gdyz_p  - gdyz_m)  / 4.d0 / eps
      dxgxz = (gdxz_p  - gdxz_m)  / 4.d0 / eps

      ax = - 0.5d0 * (gutt_p - gutt_m) / 2.d0 / eps / gutt

      bxx = ((- gutx_p / gutt_p) - (- gutx_m / gutt_m)) / 4.d0 / eps
      bxy = ((- guty_p / gutt_p) - (- guty_m / gutt_m)) / 4.d0 / eps
      bxz = ((- gutz_p / gutt_p) - (- gutz_m / gutt_m)) / 4.d0 / eps

C     Calculate y-derivatives.

      call Exact__metric(
     $     decoded_exact_model,
     $     x, y+eps, z, t,
     $     gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $     gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $     gutt_p, gutx_p, guty_p, gutz_p,
     $     guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p)
      call Exact__metric(
     $     decoded_exact_model,
     $     x, y-eps, z, t,
     $     gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $     gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $     gutt_m, gutx_m, guty_m, gutz_m,
     $     guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m)

      dygxx = (gdxx_p  - gdxx_m)  / 4.d0 / eps
      dygyy = (gdyy_p  - gdyy_m)  / 4.d0 / eps
      dygzz = (gdzz_p  - gdzz_m)  / 4.d0 / eps
      dygxy = (gdxy_p  - gdxy_m)  / 4.d0 / eps
      dygyz = (gdyz_p  - gdyz_m)  / 4.d0 / eps
      dygxz = (gdxz_p  - gdxz_m)  / 4.d0 / eps

      ay = - 0.5d0 * (gutt_p - gutt_m) / 2.d0 / eps / gutt

      byx = ((- gutx_p / gutt_p) - (- gutx_m / gutt_m)) / 4.d0 / eps
      byy = ((- guty_p / gutt_p) - (- guty_m / gutt_m)) / 4.d0 / eps
      byz = ((- gutz_p / gutt_p) - (- gutz_m / gutt_m)) / 4.d0 / eps

C     Calculate z-derivatives.

      call Exact__metric(
     $     decoded_exact_model,
     $     x, y, z+eps, t,
     $     gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $     gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $     gutt_p, gutx_p, guty_p, gutz_p,
     $     guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p)
      call Exact__metric(
     $     decoded_exact_model,
     $     x, y, z-eps, t,
     $     gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $     gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $     gutt_m, gutx_m, guty_m, gutz_m,
     $     guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m)

      dzgxx = (gdxx_p  - gdxx_m)  / 4.d0 / eps
      dzgyy = (gdyy_p  - gdyy_m)  / 4.d0 / eps
      dzgzz = (gdzz_p  - gdzz_m)  / 4.d0 / eps
      dzgxy = (gdxy_p  - gdxy_m)  / 4.d0 / eps
      dzgyz = (gdyz_p  - gdyz_m)  / 4.d0 / eps
      dzgxz = (gdxz_p  - gdxz_m)  / 4.d0 / eps

      az = - 0.5d0 * (gutt_p - gutt_m) / 2.d0 / eps / gutt

      bzx = ((- gutx_p / gutt_p) - (- gutx_m / gutt_m)) / 4.d0 / eps
      bzy = ((- guty_p / gutt_p) - (- guty_m / gutt_m)) / 4.d0 / eps
      bzz = ((- gutz_p / gutt_p) - (- gutz_m / gutt_m)) / 4.d0 / eps

C     Calculate t-derivatives, and extrinsic curvature.

      call Exact__metric(
     $     decoded_exact_model,
     $     x, y, z, t+eps,
     $     gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $     gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $     gutt_p, gutx_p, guty_p, gutz_p,
     $     guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p)
      call Exact__metric(
     $     decoded_exact_model,
     $     x, y, z, t-eps,
     $     gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $     gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $     gutt_m, gutx_m, guty_m, gutz_m,
     $     guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m)

      hxx = - (gdxx_p - gdxx_m) / 4.d0 / eps / alp
     $     + (dxgxx * betax + dygxx * betay + dzgxx * betaz
     $     + 2.d0 * (bxx * gxx + bxy * gxy + bxz * gxz)) / alp

      hyy = - (gdyy_p - gdyy_m) / 4.d0 / eps / alp
     $     + (dxgyy * betax + dygyy * betay + dzgyy * betaz
     $     + 2.d0 * (byx * gxy + byy * gyy + byz * gyz)) / alp

      hzz = - (gdzz_p - gdzz_m) / 4.d0 / eps / alp
     $     + (dxgzz * betax + dygzz * betay + dzgzz * betaz
     $     + 2.d0 * (bzx * gxz + bzy * gyz + bzz * gzz)) / alp

      hxy = - (gdxy_p - gdxy_m) / 4.d0 / eps / alp
     $     + (dxgxy * betax + dygxy * betay + dzgxy * betaz
     $     + bxx * gxy + bxy * gyy + bxz * gyz
     $     + byx * gxx + byy * gxy + byz * gxz) / alp

      hyz = - (gdyz_p - gdyz_m) / 4.d0 / eps / alp
     $     + (dxgyz * betax + dygyz * betay + dzgyz * betaz
     $     + byx * gxz + byy * gyz + byz * gzz
     $     + bzx * gxy + bzy * gyy + bzz * gyz) / alp

      hxz = - (gdxz_p - gdxz_m) / 4.d0 / eps / alp
     $     + (dxgxz * betax + dygxz * betay + dzgxz * betaz      
     $     + bxx * gxz + bxy * gyz + bxz * gzz
     $     + bzx * gxx + bzy * gxy + bzz * gxz) / alp

      return
      end
      

