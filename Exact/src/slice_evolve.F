C     Evolve the slice in the exact spacetime.

#include "cctk.h"
#include "cctk_arguments.h"

      subroutine slice_evolve(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS

      integer i,j,k,m,n
      integer nx,ny,nz

      CCTK_REAL s1d(4,3), nd(4), nu(4), norm, gd(4,4), gu(4,4)
      CCTK_REAL dx,dy,dz,dt

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      dx = cctk_delta_space(1)
      dy = cctk_delta_space(2)
      dz = cctk_delta_space(3)

      dt = cctk_delta_time

C     Lax, or forward in time, step. 
C     dxA/dt has previously been stored in slicetmp2
C     by slice_data.

      slicetmp1x = slicex + 0.5d0 * dt * slicetmp2x
      slicetmp1y = slicey + 0.5d0 * dt * slicetmp2y
      slicetmp1z = slicez + 0.5d0 * dt * slicetmp2z
      slicetmp1t = slicet + 0.5d0 * dt * slicetmp2t

C     Synchronize and bound slice.

      call linextraponebound(CCTK_FARGUMENTS,slicetmp1x)
      call synconefunc(slicetmp1x)
      call linextraponebound(CCTK_FARGUMENTS,slicetmp1y)
      call synconefunc(slicetmp1y)
      call linextraponebound(CCTK_FARGUMENTS,slicetmp1z)
      call synconefunc(slicetmp1z)
      call linextraponebound(CCTK_FARGUMENTS,slicetmp1t)
      call synconefunc(slicetmp1t)

C     Prepare leapfrog step.
C     Sum over interior points on the slice, now at midpoint slicetmp1.

      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1

C     Calculate first derivatives of slice coordinates.

               s1d(1,1) = 0.5d0*(slicetmp1x(i+1,j,k) - slicetmp1x(i-1,j,k))/dx
               s1d(1,2) = 0.5d0*(slicetmp1x(i,j+1,k) - slicetmp1x(i,j-1,k))/dy
               s1d(1,3) = 0.5d0*(slicetmp1x(i,j,k+1) - slicetmp1x(i,j,k-1))/dz

               s1d(2,1) = 0.5d0*(slicetmp1y(i+1,j,k) - slicetmp1y(i-1,j,k))/dx
               s1d(2,2) = 0.5d0*(slicetmp1y(i,j+1,k) - slicetmp1y(i,j-1,k))/dy
               s1d(2,3) = 0.5d0*(slicetmp1y(i,j,k+1) - slicetmp1y(i,j,k-1))/dz

               s1d(3,1) = 0.5d0*(slicetmp1z(i+1,j,k) - slicetmp1z(i-1,j,k))/dx
               s1d(3,2) = 0.5d0*(slicetmp1z(i,j+1,k) - slicetmp1z(i,j-1,k))/dy
               s1d(3,3) = 0.5d0*(slicetmp1z(i,j,k+1) - slicetmp1z(i,j,k-1))/dz

               s1d(4,1) = 0.5d0*(slicetmp1t(i+1,j,k) - slicetmp1t(i-1,j,k))/dx
               s1d(4,2) = 0.5d0*(slicetmp1t(i,j+1,k) - slicetmp1t(i,j-1,k))/dy
               s1d(4,3) = 0.5d0*(slicetmp1t(i,j,k+1) - slicetmp1t(i,j,k-1))/dz

C     Now we need the exact solution metric in the preferred coordinates
C     x^A. 
               call exactmetric(
     $              slicetmp1x(i,j,k), slicetmp1y(i,j,k), slicetmp1z(i,j,k), 
     $              slicetmp1t(i,j,k),
     $              gd(4,4), gd(1,4), gd(2,4), gd(3,4),
     $              gd(1,1), gd(2,2), gd(3,3), 
     $              gd(1,2), gd(2,3), gd(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3))

C     Calculate n^A and dx^A/dt
#include "include/slice_normal.h"

            end do
         end do
      end do

C     Synchronize and bound slicetmp2, which contains dx^A/dt.

      call linextraponebound(CCTK_FARGUMENTS,slicetmp2x)
      call synconefunc(slicetmp2x)
      call linextraponebound(CCTK_FARGUMENTS,slicetmp2y)
      call synconefunc(slicetmp2y)
      call linextraponebound(CCTK_FARGUMENTS,slicetmp2z)
      call synconefunc(slicetmp2z)
      call linextraponebound(CCTK_FARGUMENTS,slicetmp2t)
      call synconefunc(slicetmp2t)
 
C Leapfrog step.

      slicex = slicex + dt * slicetmp2x
      slicey = slicey + dt * slicetmp2y
      slicez = slicez + dt * slicetmp2z
      slicet = slicet + dt * slicetmp2t

C     Synchronize and bound slice.

      call linextraponebound(CCTK_FARGUMENTS,slicex)
      call synconefunc(slicex)
      call linextraponebound(CCTK_FARGUMENTS,slicey)
      call synconefunc(slicey)
      call linextraponebound(CCTK_FARGUMENTS,slicez)
      call synconefunc(slicez)
      call linextraponebound(CCTK_FARGUMENTS,slicet)
      call synconefunc(slicet)

C     Extract Cauchy data at the new position, and store dxA/dt 
C     for use in the next Lax step.

      call slice_data(CCTK_FARGUMENTS)

      return
      end



