#include "cctk.h"
#include "cctk_parameters.h"
#include "cctk_arguments.h"

      subroutine exactgauge(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k
      integer nx,ny,nz
      integer CCTK_Equals

      CCTK_REAL tplushalf
      CCTK_REAL gxxjunk, gyyjunk, gzzjunk, 
     $     gxyjunk, gyzjunk, gxzjunk,
     $     hxxjunk, hyyjunk, hzzjunk, 
     $     hxyjunk, hyzjunk, hxzjunk,
     $     dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $     dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $     dygxxjunk, dygyyjunk, dygzzjunk, 
     $     dygxyjunk, dygyzjunk, dygxzjunk,
     $     dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $     dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $     alpjunk, axjunk, ayjunk, azjunk, 
     $     betaxjunk, betayjunk, betazjunk,
     $     bxxjunk, bxyjunk, bxzjunk, 
     $     byxjunk, byyjunk, byzjunk, 
     $     bzxjunk, bzyjunk, bzzjunk

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     Set t = time + dt/2.  This is for consistency with a second
C     order numerical scheme.

      tplushalf = cctk_time + 0.5D0*cctk_delta_time

C     Set both lapse and shift.

      if ((CCTK_Equals(slicing,"exact").ne.0).and.
     .    (CCTK_Equals(shift,"exact").ne.0)) then

         do k=1,nz
            do j=1,ny
               do i=1,nx

                  call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), tplushalf,
     $                 gxxjunk, gyyjunk, gzzjunk, 
     $                 gxyjunk, gyzjunk, gxzjunk,
     $                 hxxjunk, hyyjunk, hzzjunk, 
     $                 hxyjunk, hyzjunk, hxzjunk,
     $                 dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $                 dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $                 dygxxjunk, dygyyjunk, dygzzjunk, 
     $                 dygxyjunk, dygyzjunk, dygxzjunk,
     $                 dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $                 dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $                 alp(i,j,k), axjunk, ayjunk, azjunk, 
     $                 betax(i,j,k), betay(i,j,k), betaz(i,j,k),
     $                 bxxjunk, bxyjunk, bxzjunk, 
     $                 byxjunk, byyjunk, byzjunk, 
     $                 bzxjunk, bzyjunk, bzzjunk)

               end do
            end do
         end do

C     Set lapse only.

      else if (CCTK_Equals(slicing,"exact").ne.0) then

         do k=1,nz
            do j=1,ny
               do i=1,nx

                  call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), tplushalf,
     $                 gxxjunk, gyyjunk, gzzjunk, 
     $                 gxyjunk, gyzjunk, gxzjunk,
     $                 hxxjunk, hyyjunk, hzzjunk, 
     $                 hxyjunk, hyzjunk, hxzjunk,
     $                 dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $                 dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $                 dygxxjunk, dygyyjunk, dygzzjunk, 
     $                 dygxyjunk, dygyzjunk, dygxzjunk,
     $                 dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $                 dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $                 alp(i,j,k), axjunk, ayjunk, azjunk, 
     $                 betaxjunk, betayjunk, betazjunk,
     $                 bxxjunk, bxyjunk, bxzjunk, 
     $                 byxjunk, byyjunk, byzjunk, 
     $                 bzxjunk, bzyjunk, bzzjunk)

               end do
            end do
         end do

C     Set shift only.

      else if (CCTK_Equals(shift,"exact").ne.0) then

         do k=1,nz
            do j=1,ny
               do i=1,nx

                  call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), tplushalf,
     $                 gxxjunk, gyyjunk, gzzjunk, 
     $                 gxyjunk, gyzjunk, gxzjunk,
     $                 hxxjunk, hyyjunk, hzzjunk, 
     $                 hxyjunk, hyzjunk, hxzjunk,
     $                 dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $                 dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $                 dygxxjunk, dygyyjunk, dygzzjunk, 
     $                 dygxyjunk, dygyzjunk, dygxzjunk,
     $                 dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $                 dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $                 alpjunk, axjunk, ayjunk, azjunk, 
     $                 betax(i,j,k), betay(i,j,k), betaz(i,j,k),
     $                 bxxjunk, bxyjunk, bxzjunk, 
     $                 byxjunk, byyjunk, byzjunk, 
     $                 bzxjunk, bzyjunk, bzzjunk)

               end do
            end do
         end do

      end if

      return
      end
