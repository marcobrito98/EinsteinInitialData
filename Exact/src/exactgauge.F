C     This routine sets the lapse and/or shift by calling a routine
C     that does it pointwise. Note that it could be easily modified
C     to set the Bona-Masso variables B_xx etc. 

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine exactgauge(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
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

      if (((CCTK_Equals(lapse_evolution_method,"exact").ne.0).and.
     $     (CCTK_Equals(shift_evolution_method,"exact").ne.0))
     $     .or.
     $     ((CCTK_Equals(initial_lapse,"exact").ne.0).and.
     $     (CCTK_Equals(initial_shift,"exact").ne.0))) then

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

      else if ((CCTK_Equals(lapse_evolution_method,"exact").ne.0)
     $        .or. (CCTK_Equals(initial_lapse,"exact").ne.0)) then

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

      else if ((CCTK_Equals(shift_evolution_method,"exact").ne.0) 
     $        .or. (CCTK_Equals(initial_shift,"exact").ne.0)) then

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

      else
         call CCTK_WARN(1,'exactgauge has been called without doing anything')
      end if

      return
      end
