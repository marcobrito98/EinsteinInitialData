
#include "cctk.h"
#include "cctk_parameters.h"
#include "cctk_arguments.h"

      subroutine exactboundary(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS

      integer i,j,k
      integer nx,ny,nz

      CCTK_REAL tplusone
      CCTK_REAL 
     $     dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $     dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $     dygxxjunk, dygyyjunk, dygzzjunk, 
     $     dygxyjunk, dygyzjunk, dygxzjunk,
     $     dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $     dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $     axjunk, ayjunk, azjunk, 
     $     bxxjunk, bxyjunk, bxzjunk, 
     $     byxjunk, byyjunk, byzjunk, 
     $     bzxjunk, bzyjunk, bzzjunk

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     Set all initial data including dijk and vi on all points which 
C     are on the boundary of the domain if it really is the boundary
C     of the complete grid. Treat all six sides of the grid cube this way.

c     Set t = time + dt.  This is necessary here because by the time
c     we reach this point the geometry has been evolved one time step
c     but the variable `time' still hasn't been updated. 

      tplusone = cctk_time + cctk_delta_time

C     Note we also always set the lapse and shift at the boundaries at
C     time t+1. This is to provide boundary conditions for testing
C     elliptic gauge conditions. If they are not used, they will be
C     overwritten by exactgauge.

#define EXACTDATAPOINT \
               call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), tplusone, \
                   gxx(i,j,k), gyy(i,j,k), gzz(i,j,k), \
                   gxy(i,j,k), gyz(i,j,k), gxz(i,j,k), \
                   kxx(i,j,k), kyy(i,j,k), kzz(i,j,k), \
                   kxy(i,j,k), kyz(i,j,k), kxz(i,j,k), \
                   dxgxxjunk, dxgyyjunk, dxgzzjunk, \
                   dxgxyjunk, dxgyzjunk, dxgxzjunk, \
                   dygxxjunk, dygyyjunk, dygzzjunk, \
                   dygxyjunk, dygyzjunk, dygxzjunk, \
                   dzgxxjunk, dzgyyjunk, dzgzzjunk, \
                   dzgxyjunk, dzgyzjunk, dzgxzjunk, \
                   alp(i,j,k), axjunk, ayjunk, azjunk, \
                   betax(i,j,k), betay(i,j,k), betaz(i,j,k), \
                   bxxjunk, bxyjunk, bxzjunk, \
                   byxjunk, byyjunk, byzjunk, \
                   bzxjunk, bzyjunk, bzzjunk) &&\

      if (cctk_bbox(1) .eq. 1) then 
         i=1
         do j=1,ny
            do k=1,nz
               EXACTDATAPOINT
            end do
         end do
      end if

      if (cctk_bbox(2) .eq. 1) then 
         i=nx
         do j=1,ny
            do k=1,nz
               EXACTDATAPOINT
            end do
         end do
      end if

      if (cctk_bbox(3) .eq. 1) then 
         j=1
         do i=1,nx
            do k=1,nz
               EXACTDATAPOINT
            end do
         end do
      end if

      if (cctk_bbox(4) .eq. 1) then 
         j=ny
         do i=1,nx
            do k=1,nz
               EXACTDATAPOINT
           end do
         end do
      end if

      if (cctk_bbox(5) .eq. 1) then 
         k=1
         do j=1,ny
            do i=1,nx
               EXACTDATAPOINT
            end do
         end do
      end if

      if (cctk_bbox(6) .eq. 1) then 
         k=nz
         do j=1,ny
            do i=1,nx
               EXACTDATAPOINT
            end do
         end do
      end if

      return
      end

