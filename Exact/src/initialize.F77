C     Wrapper for boostrotdata. Calls it and vectorini.
C     Sets Cauchy data, lapse and shift, and what else is needed
C     in the Bona-Masso formalism, at an initial time.
C $Header$
      
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine Exact__initialize(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer nx,ny,nz

      CCTK_REAL tt, xx, yy, zz
      CCTK_REAL alpjunk, dtalpjunk, axjunk, ayjunk, azjunk, 
     $     betaxjunk, betayjunk, betazjunk,
     $     dtbetaxjunk, dtbetayjunk, dtbetazjunk,
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
      CCTK_REAL
     $     exact_psi,
     $     exact_psix, exact_psiy, exact_psiz,
     $     exact_psixx, exact_psiyy, exact_psizz,
     $     exact_psixy, exact_psiyz, exact_psixz

      call CCTK_INFO('setting exact data on slice')

C     Set conformal state
      if (CCTK_EQUALS(metric_type, "static conformal")) then
        conformal_state=1
        if (CCTK_EQUALS(conformal_storage,"factor+derivs")) then
          conformal_state = 2
        else if
     $     (CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs"))
     $     then
          conformal_state = 3
        end if
      end if

C     Note I assume time has been initialized to physical time.
C     Set data pointwise.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C$o mp parallel do private(
C$o mp$     i, j, k,
C$o mp$     tt, xx, yy, zz,
C$o mp$     alpjunk, dtalpjunk, axjunk, ayjunk, azjunk, 
C$o mp$     betaxjunk, betayjunk, betazjunk,
C$o mp$     dtbetaxjunk, dtbetayjunk, dtbetazjunk,
C$o mp$     bxxjunk, bxyjunk, bxzjunk, 
C$o mp$     byxjunk, byyjunk, byzjunk, 
C$o mp$     bzxjunk, bzyjunk, bzzjunk,
C$o mp$     dxgxxjunk, dxgyyjunk, dxgzzjunk, 
C$o mp$     dxgxyjunk, dxgyzjunk, dxgxzjunk,
C$o mp$     dygxxjunk, dygyyjunk, dygzzjunk, 
C$o mp$     dygxyjunk, dygyzjunk, dygxzjunk,
C$o mp$     dzgxxjunk, dzgyyjunk, dzgzzjunk, 
C$o mp$     dzgxyjunk, dzgyzjunk, dzgxzjunk,
C$o mp$     exact_psi,
C$o mp$     exact_psix, exact_psiy, exact_psiz,
C$o mp$     exact_psixx, exact_psiyy, exact_psizz,
C$o mp$     exact_psixy, exact_psiyz, exact_psixz)
      do k=1,nz
         do j=1,ny
            do i=1,nx

               tt = cctk_time
               xx = x(i,j,k) - cctk_time * shift_add_x
               yy = y(i,j,k) - cctk_time * shift_add_y
               zz = z(i,j,k) - cctk_time * shift_add_z

C     Initialize the psi of exact
C     (also to tell the models about the conformal_state)
               if (conformal_state .ne. 0) then
                 exact_psi = 1.0D0
               else
                 exact_psi = 0.0D0
               end if
               exact_psix = 0.0D0
               exact_psiy = 0.0D0
               exact_psiz = 0.0D0
               exact_psixx = 0.0D0
               exact_psiyy = 0.0D0
               exact_psizz = 0.0D0
               exact_psixy = 0.0D0
               exact_psiyz = 0.0D0
               exact_psixz = 0.0D0

               call Exact__Bona_Masso_data(
     $              decoded_exact_model,
     $              xx, yy, zz, tt,
     $              gxx(i,j,k), gyy(i,j,k), gzz(i,j,k), 
     $              gxy(i,j,k), gyz(i,j,k), gxz(i,j,k),
     $              kxx(i,j,k), kyy(i,j,k), kzz(i,j,k), 
     $              kxy(i,j,k), kyz(i,j,k), kxz(i,j,k),
     $              exact_psi,
     $              exact_psix, exact_psiy, exact_psiz,
     $              exact_psixx, exact_psiyy, exact_psizz,
     $              exact_psixy, exact_psiyz, exact_psixz,
     $              dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $              dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $              dygxxjunk, dygyyjunk, dygzzjunk, 
     $              dygxyjunk, dygyzjunk, dygxzjunk,
     $              dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $              dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $              alpjunk, dtalpjunk, axjunk, ayjunk, azjunk, 
     $              betaxjunk, betayjunk, betazjunk,
     $              dtbetaxjunk, dtbetayjunk, dtbetazjunk,
     $              bxxjunk, bxyjunk, bxzjunk,
     $              byxjunk, byyjunk, byzjunk,
     $              bzxjunk, bzyjunk, bzzjunk)

C     Save the conformal factor if wanted
               if (conformal_state .ne. 0) then
                 psi(i,j,k) = exact_psi
                 if (conformal_state .gt. 1) then
                   psix(i,j,k) = exact_psix
                   psiy(i,j,k) = exact_psiy
                   psiz(i,j,k) = exact_psiz
                   if (conformal_state .gt. 2) then
                     psixx(i,j,k) = exact_psixx
                     psiyy(i,j,k) = exact_psiyy
                     psizz(i,j,k) = exact_psizz
                     psixy(i,j,k) = exact_psixy
                     psiyz(i,j,k) = exact_psiyz
                     psixz(i,j,k) = exact_psixz
                   end if
                 end if
               end if

            end do
         end do
      end do
      
C     Tell the code there is no need to treat the conformal factor
C     as a separate field. That is, we have set the physical metric here.
c     Commented out in einstein revamp, now Exact does not inherit anything
c     about the conformal factor
c     Now it does again (see above, knarf)

      return
      end
