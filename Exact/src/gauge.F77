C     This routine sets the lapse and/or shift by calling a routine
C     that does it pointwise. Note that it could be easily modified
C     to set the Bona-Masso variables B_xx etc. 
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine Exact__gauge(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer nx,ny,nz

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
      LOGICAL is_initial_slice, is_later_slice

C     are we on the initial slice or some later slice?
C     n.b. the logical expressions later in this function involving
C          these flags below would be *so* much nicer if Fortran
C          grokked C-style conditional expressions... :) :)
      is_initial_slice = cctk_iteration .eq. 0
      is_later_slice   = cctk_iteration .ne. 0

C     Grid parameters.
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     This code used to set  t = time + dt/2  to get 2nd order accuracy,
C     but this leads to the initial data being set at the wrong time. :(
C     In the context of MoL, we want to set variables at the standard Cactus
C     time (cctk_time), because MoL takes care of calling us at each MoL
C     iteration, and updating the field variables appropriately.
C
C     Alas, setting at cctk_time probably gives O(dt) errors for non-MoL
C     evoutions where Exact is used to set stuff at each time step.
C     Fixing this [unless we just declare all non-MoL stuff obselete :) ]
C     probably requires cleaning up our (++messy)  schedule.ccl , which
C     is why this remains a bug for now... :( :(

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C     Set both lapse and shift?
C
      if (
     $     ( is_initial_slice
     $       .and. (CCTK_Equals(initial_lapse,         "exact").ne.0)
     $       .and. (CCTK_Equals(initial_shift,         "exact").ne.0) )
     $     .or.
     $     ( is_later_slice
     $       .and. (CCTK_Equals(lapse_evolution_method,"exact").ne.0)
     $       .and. (CCTK_Equals(shift_evolution_method,"exact").ne.0) )
     $   ) then

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
     $                 x(i,j,k), y(i,j,k), z(i,j,k), cctk_time,
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
                  betax(i,j,k) = betax(i,j,k) + shift_add_x
                  betay(i,j,k) = betay(i,j,k) + shift_add_y
                  betaz(i,j,k) = betaz(i,j,k) + shift_add_z
               end do
            end do
         end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C     Set lapse only?
C
      elseif (
     $         ( is_initial_slice
     $           .and. (CCTK_Equals(initial_lapse,         "exact").ne.0) )
     $         .or.
     $         ( is_later_slice
     $           .and. (CCTK_Equals(lapse_evolution_method,"exact").ne.0) )
     $       ) then

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
     $                 x(i,j,k), y(i,j,k), z(i,j,k), cctk_time,
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C     Set shift only?
C
      elseif (
     $         ( is_initial_slice
     $           .and. (CCTK_Equals(initial_shift,"exact").ne.0) )
     $         .or.
     $         ( is_later_slice
     $           .and. (CCTK_Equals(shift_evolution_method,"exact").ne.0) )
     $       ) then

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
     $                 x(i,j,k), y(i,j,k), z(i,j,k), cctk_time,
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
                  betax(i,j,k) = betax(i,j,k) + shift_add_x
                  betay(i,j,k) = betay(i,j,k) + shift_add_y
                  betaz(i,j,k) = betaz(i,j,k) + shift_add_z

               end do
            end do
         end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      else
         call CCTK_WARN(1,'Exact__gauge has been called without doing anything')
      end if

      return
      end
