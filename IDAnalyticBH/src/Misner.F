C     This wrapper written by Carsten Gundlach. 

#include "cctk.h"
#include "declare_arguments.h"
#include "declare_parameters.h"

      subroutine Misner(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_PARAMETERS

      integer i, j, k
      REAL one, zero, tmp0, tmp1, tmp2, tmp3, tmp4, nm_eps

C     Finite differencing step.
      parameter(nm_eps = 1.d-6)

      one = 1.0D0
      zero = 0.0D0

c     Must use conformal metric
c     -------------------------
      if (conformal_state == 0) then
         call CCTK_Warn(1,"Changing to conformal metric in Misner")
         conformal_state = 1
      end if

C     Initialize some C global variables before calls to 
C     misner_psi_eval.
      call nmisner_init(misner_nbh,mu,nmax)

C     Get value of psi pointwise from a C function.
      do k=1,sh(1)
         do j=1,sh(2)
            do i=1,sh(3)
 
               call nmisner_eval_psi(x(i,j,k),y(i,j,k),z(i,j,k),tmp0)
               psi(i,j,k) = tmp0
               
               call nmisner_eval_psi(x(i,j,k)+nm_eps,y(i,j,k),z(i,j,k),tmp1)
               call nmisner_eval_psi(x(i,j,k)-nm_eps,y(i,j,k),z(i,j,k),tmp2)
               psix(i,j,k) = 0.5d0 * (tmp1 - tmp2) / nm_eps
               psixx(i,j,k) = (tmp1 + tmp2 - 2.d0 * tmp0) / nm_eps**2

               call nmisner_eval_psi(x(i,j,k),y(i,j,k)+nm_eps,z(i,j,k),tmp1)
               call nmisner_eval_psi(x(i,j,k),y(i,j,k)-nm_eps,z(i,j,k),tmp2)
               psiy(i,j,k) = 0.5d0 * (tmp1 - tmp2) / nm_eps
               psiyy(i,j,k) = (tmp1 + tmp2 - 2.d0 * tmp0) / nm_eps**2

               call nmisner_eval_psi(x(i,j,k),y(i,j,k),z(i,j,k)+nm_eps,tmp1)
               call nmisner_eval_psi(x(i,j,k),y(i,j,k),z(i,j,k)-nm_eps,tmp2)
               psiz(i,j,k) = 0.5d0 * (tmp1 - tmp2) / nm_eps
               psizz(i,j,k) = (tmp1 + tmp2 - 2.d0 * tmp0) / nm_eps**2

               call nmisner_eval_psi(x(i,j,k)+nm_eps,y(i,j,k)+nm_eps,z(i,j,k),tmp1)
               call nmisner_eval_psi(x(i,j,k)+nm_eps,y(i,j,k)-nm_eps,z(i,j,k),tmp2)
               call nmisner_eval_psi(x(i,j,k)-nm_eps,y(i,j,k)+nm_eps,z(i,j,k),tmp3)
               call nmisner_eval_psi(x(i,j,k)-nm_eps,y(i,j,k)-nm_eps,z(i,j,k),tmp4)
               psixy(i,j,k) = 0.25d0 * (tmp1 - tmp2 - tmp3 + tmp4) / nm_eps**2

               call nmisner_eval_psi(x(i,j,k),y(i,j,k)+nm_eps,z(i,j,k)+nm_eps,tmp1)
               call nmisner_eval_psi(x(i,j,k),y(i,j,k)-nm_eps,z(i,j,k)+nm_eps,tmp2)
               call nmisner_eval_psi(x(i,j,k),y(i,j,k)+nm_eps,z(i,j,k)-nm_eps,tmp3)
               call nmisner_eval_psi(x(i,j,k),y(i,j,k)-nm_eps,z(i,j,k)-nm_eps,tmp4)
               psiyz(i,j,k) = 0.25d0 * (tmp1 - tmp2 - tmp3 + tmp4) / nm_eps**2

               call nmisner_eval_psi(x(i,j,k)+nm_eps,y(i,j,k),z(i,j,k)+nm_eps,tmp1)
               call nmisner_eval_psi(x(i,j,k)+nm_eps,y(i,j,k),z(i,j,k)-nm_eps,tmp2)
               call nmisner_eval_psi(x(i,j,k)-nm_eps,y(i,j,k),z(i,j,k)+nm_eps,tmp3)
               call nmisner_eval_psi(x(i,j,k)-nm_eps,y(i,j,k),z(i,j,k)-nm_eps,tmp4)
               psixz(i,j,k) = 0.25d0 * (tmp1 - tmp2 - tmp3 + tmp4) / nm_eps**2
            end do
         end do
      end do

C     Cactus conventions:
      psix = psix / psi
      psiy = psiy / psi
      psiz = psiz / psi
      psixx = psixx / psi
      psixy = psixy / psi
      psixz = psixz / psi
      psiyy = psiyy / psi
      psiyz = psiyz / psi
      psizz = psizz / psi

C     The conformal metric is flat.
      gxx = one
      gyy = one
      gzz = one
      gxy = zero 
      gxz = zero  
      gyz = zero   

C     Time-symmetric data.
      hxx = zero
      hyy = zero
      hzz = zero
      hxy = zero 
      hxz = zero  
      hyz = zero   

      return
      end

