C     This routine sets the lapse and/or shift by calling a routine
C     that does it pointwise. Note that it could be easily modified
C     to set the Bona-Masso variables B_xx etc. 
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine Exact__gauge(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer nx,ny,nz

      CCTK_REAL t_exact
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
      CCTK_REAL
     $     exact_psi,
     $     exact_psix, exact_psiy, exact_psiz,
     $     exact_psixx, exact_psiyy, exact_psizz,
     $     exact_psixy, exact_psiyz, exact_psixz

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     This code used to set  t = time + dt/2  to get 2nd order accuracy,
C     but this leads to the initial data being set at the wrong time. :(
C     In the context of MoL, we want to set the "exact" time anyway, because
C     MoL takes care of calling us at each MoL iteration, and updating the
C     field variables appropriately.
C
C     Alas, the current version of this code probably gives O(dt) errors
C     for non-MoL evoutions where Exact is used to set stuff at each time
C     step.  Fixing this probably requires cleaning up our (++messy)
C      schedule.ccl , which is why this remains a bug for now... :( :(

      t_exact = cctk_time

C     Set both lapse and shift.

      if (( (CCTK_Equals(lapse_evolution_method,"exact").ne.0) .and.
     $      (CCTK_Equals(shift_evolution_method,"exact").ne.0) )
     $    .or.
     $    ( (CCTK_Equals(initial_lapse,"exact").ne.0) .and.
     $      (CCTK_Equals(initial_shift,"exact").ne.0) )) then

         do k=1,nz
            do j=1,ny
               do i=1,nx

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
                  exact_psixx = 0.0D0
                  exact_psizz = 0.0D0
                  exact_psixy = 0.0D0
                  exact_psiyz = 0.0D0
                  exact_psixz = 0.0D0

                  call Exact__Bona_Masso_data(
     $                 decoded_exact_model,
     $                 x(i,j,k), y(i,j,k), z(i,j,k), t_exact,
     $                 gxxjunk, gyyjunk, gzzjunk, 
     $                 gxyjunk, gyzjunk, gxzjunk,
     $                 hxxjunk, hyyjunk, hzzjunk, 
     $                 hxyjunk, hyzjunk, hxzjunk,
     $                 exact_psi,
     $                 exact_psix, exact_psiy, exact_psiz,
     $                 exact_psixx, exact_psiyy, exact_psizz,
     $                 exact_psixy, exact_psiyz, exact_psixz,
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
                  exact_psixx = 0.0D0
                  exact_psizz = 0.0D0
                  exact_psixy = 0.0D0
                  exact_psiyz = 0.0D0
                  exact_psixz = 0.0D0

                  call Exact__Bona_Masso_data(
     $                 decoded_exact_model,
     $                 x(i,j,k), y(i,j,k), z(i,j,k), t_exact,
     $                 gxxjunk, gyyjunk, gzzjunk, 
     $                 gxyjunk, gyzjunk, gxzjunk,
     $                 hxxjunk, hyyjunk, hzzjunk, 
     $                 hxyjunk, hyzjunk, hxzjunk,
     $                 exact_psi,
     $                 exact_psix, exact_psiy, exact_psiz,
     $                 exact_psixx, exact_psiyy, exact_psizz,
     $                 exact_psixy, exact_psiyz, exact_psixz,
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
                  exact_psixx = 0.0D0
                  exact_psizz = 0.0D0
                  exact_psixy = 0.0D0
                  exact_psiyz = 0.0D0
                  exact_psixz = 0.0D0

                  call Exact__Bona_Masso_data(
     $                 decoded_exact_model,
     $                 x(i,j,k), y(i,j,k), z(i,j,k), t_exact,
     $                 gxxjunk, gyyjunk, gzzjunk, 
     $                 gxyjunk, gyzjunk, gxzjunk,
     $                 hxxjunk, hyyjunk, hzzjunk, 
     $                 hxyjunk, hyzjunk, hxzjunk,
     $                 exact_psi,
     $                 exact_psix, exact_psiy, exact_psiz,
     $                 exact_psixx, exact_psiyy, exact_psizz,
     $                 exact_psixy, exact_psiyz, exact_psixz,
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
         call CCTK_WARN(1,'Exact__gauge has been called without doing anything')
      end if

      return
      end
