C     Wrapper for boostrotdata. Calls it and vectorini.
C     Sets Cauchy data, lapse and shift, and what else is needed
C     in the Bona-Masso formalism, at an initial time.
      
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine exactinitialize(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k
      integer nx,ny,nz
      integer CCTK_Equals

      CCTK_REAL alpjunk, axjunk, ayjunk, azjunk, 
     $     betaxjunk, betayjunk, betazjunk,
     $     bxxjunk, bxyjunk, bxzjunk, 
     $     byxjunk, byyjunk, byzjunk, 
     $     bzxjunk, bzyjunk, bzzjunk
      CCTK_REAL 
     $     dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $     dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $     dygxxjunk, dygyyjunk, dygzzjunk, 
     $     dygxyjunk, dygyzjunk, dygxzjunk,
     $     dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $     dzgxyjunk, dzgyzjunk, dzgxzjunk

C     Note I assume time has been initialized to physical time.
C     Set data pointwise.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      do k=1,nz
         do j=1,ny
            do i=1,nx
               
               call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), cctk_time,
     $              gxx(i,j,k), gyy(i,j,k), gzz(i,j,k), 
     $              gxy(i,j,k), gyz(i,j,k), gxz(i,j,k),
     $              kxx(i,j,k), kyy(i,j,k), kzz(i,j,k), 
     $              kxy(i,j,k), kyz(i,j,k), kxz(i,j,k),
     $              dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $              dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $              dygxxjunk, dygyyjunk, dygzzjunk, 
     $              dygxyjunk, dygyzjunk, dygxzjunk,
     $              dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $              dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $              alpjunk, axjunk, ayjunk, azjunk, 
     $              betaxjunk, betayjunk, betazjunk,
     $              bxxjunk, bxyjunk, bxzjunk,
     $              byxjunk, byyjunk, byzjunk,
     $              bzxjunk, bzyjunk, bzzjunk)

            end do
         end do
      end do
      
C     Tell the code there is no need to treat the conformal factor
C     as a separate field. That is, we have set the physical metric here.
c     Commented out in einstein revamp, now Exact doesnot inherit anything
c     about the conformal factor
c      conformal_state = 0
c      psi = 1.0D0

      return
      end
