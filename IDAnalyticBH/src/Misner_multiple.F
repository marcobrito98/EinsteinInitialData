 /*@@
   @file      Misner_multiple.F
   @date      
   @author    Carsten Gundlach
   @desc 
      Set up initial data for multiple Misner black holes
   @enddesc 
 @@*/

#include "cctk.h"
#include "declare_arguments.h"
#include "declare_parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

 /*@@
   @routine    Misner_multiple
   @date       
   @author     Carsten Gundlach
   @desc 
      Set up initial data for multiple Misner black holes
   @enddesc 
   @calls     
   @history 
 
   @endhistory 

@@*/
      subroutine Misner_multiple(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_PARAMETERS

      integer :: i, j, k
      CCTK_REAL :: xval, yval, zval
      CCTK_REAL :: one, zero, tmp0, tmp1, tmp2, tmp3, tmp4
      CCTK_REAL, parameter :: nm_eps = 1.d-6   ! finite differencing step

      one  = 1.0D0
      zero = 0.0D0

c     Initialize C global variables
c     -----------------------------
      call Misner_init(misner_nbh,mu,nmax)

C     Get value of psi pointwise from a C function
c     --------------------------------------------
      do k=1,cctk_lsh(1)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(3)
 
               xval = x(i,j,k)
               yval = y(i,j,k)
               zval = z(i,j,k)

               call MisnerEvalPsi(xval,yval,zval,tmp0)
               psi(i,j,k) = tmp0
               
c              Only calculate derivatives of psi if required
c              ---------------------------------------------
               if (use_conformal_derivs == 1) then

                 call MisnerEvalPsi(xval+nm_eps,yval,zval,tmp1)
                 call MisnerEvalPsi(xval-nm_eps,yval,zval,tmp2)
                 psix(i,j,k) = 0.5d0*(tmp1-tmp2)/nm_eps
                 psixx(i,j,k) = (tmp1+tmp2-2.d0*tmp0)/nm_eps**2

                 call MisnerEvalPsi(xval,yval+nm_eps,zval,tmp1)
                 call MisnerEvalPsi(xval,yval-nm_eps,zval,tmp2)
                 psiy(i,j,k) = 0.5d0*(tmp1-tmp2)/nm_eps
                 psiyy(i,j,k) = (tmp1+tmp2-2.d0*tmp0)/nm_eps**2

                 call MisnerEvalPsi(xval,yval,zval+nm_eps,tmp1)
                 call MisnerEvalPsi(xval,yval,zval-nm_eps,tmp2)
                 psiz(i,j,k) = 0.5d0*(tmp1-tmp2)/nm_eps
                 psizz(i,j,k) = (tmp1+tmp2-2.d0*tmp0)/nm_eps**2

                 call MisnerEvalPsi(xval+nm_eps,yval+nm_eps,zval,tmp1)
                 call MisnerEvalPsi(xval+nm_eps,yval-nm_eps,zval,tmp2)
                 call MisnerEvalPsi(xval-nm_eps,yval+nm_eps,zval,tmp3)
                 call MisnerEvalPsi(xval-nm_eps,yval-nm_eps,zval,tmp4)
                 psixy(i,j,k) = 0.25d0*(tmp1-tmp2-tmp3+tmp4)/nm_eps**2

                 call MisnerEvalPsi(xval,yval+nm_eps,zval+nm_eps,tmp1)
                 call MisnerEvalPsi(xval,yval-nm_eps,zval+nm_eps,tmp2)
                 call MisnerEvalPsi(xval,yval+nm_eps,zval-nm_eps,tmp3)
                 call MisnerEvalPsi(xval,yval-nm_eps,zval-nm_eps,tmp4)
                 psiyz(i,j,k) = 0.25d0*(tmp1-tmp2-tmp3+tmp4)/nm_eps**2

                 call MisnerEvalPsi(xval+nm_eps,yval,zval+nm_eps,tmp1)
                 call MisnerEvalPsi(xval+nm_eps,yval,zval-nm_eps,tmp2)
                 call MisnerEvalPsi(xval-nm_eps,yval,zval+nm_eps,tmp3)
                 call MisnerEvalPsi(xval-nm_eps,yval,zval-nm_eps,tmp4)
                 psixz(i,j,k) = 0.25d0*(tmp1-tmp2-tmp3+tmp4)/nm_eps**2

              end if

            end do
         end do
      end do

c     Cactus conventions
c     ------------------
      if (use_conformal_derivs == 1) then

        psix = psix/psi
        psiy = psiy/psi
        psiz = psiz/psi
        psixx = psixx/psi
        psixy = psixy/psi
        psixz = psixz/psi
        psiyy = psiyy/psi
        psiyz = psiyz/psi
        psizz = psizz/psi

      end if

c     Metric depends on conformal state
c     ---------------------------------
      if (conformal_state == CONFORMAL_METRIC) then

        gxx = one
        gyy = one
        gzz = one
        gxy = zero 
        gxz = zero  
        gyz = zero   

      else 

        gxx = psi**4
        gyy = gxx
        gzz = gxx
        gxy = zero
        gxz = zero
        gyz = zero

      end if

c     Time-symmetric data
c     -------------------
      kxx = zero
      kyy = zero
      kzz = zero
      kxy = zero 
      kxz = zero  
      kyz = zero   

      return
      end

