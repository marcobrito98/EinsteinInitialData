#include "cctk.h"
#include "declare_arguments.h"
#include "declare_parameters.h"

c Need include file from Einstein
#include "../../packages/CactusEinstein/Einstein/src/Einstein.h"

      subroutine Schwarzschild(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_PARAMETERS

      REAL zero,one,two,three

      REAL t1,r0

      integer CCTK_Equals
      integer i,j,k
      integer is,js,ks,ie,je,ke
      REAL inval

      zero = 0.d0
      one = 1.d0
      two = 2.d0
      three = 3.d0
            
c     conformal metric flag
      if (use_conformal == 1) then

         conformal_state = CONFORMAL_METRIC

c        Compute conformal factor
         psi = ( one + mass/two/r)
         
c        derivatives of psi / psi
         psix = -x*mass/two/r**3 / psi
         psiy = -y*mass/two/r**3 / psi
         psiz = -z*mass/two/r**3 / psi
         
         psixy = three*x*y*mass/two/r**5 /psi
         psixz = three*x*z*mass/two/r**5 /psi
         psiyz = three*y*z*mass/two/r**5 /psi
         
         psixx  = (three*x**2 - r**2)*mass/two/r**5/psi
         psiyy  = (three*y**2 - r**2)*mass/two/r**5/psi
         psizz  = (three*z**2 - r**2)*mass/two/r**5/psi          
         
         gxx = one
         gyy = one
         gzz = one
         gxy = zero 
         gxz = zero  
         gyz = zero   
          
      else

         conformal_state = NOCONFORMAL_METRIC

         gxx = ( one + mass/two/r)**4
         gyy = gxx
         gzz = gxx
         gxy = zero
         gxz = zero
         gyz = zero
         
      endif

c     If the initial lapse is not one ...
      if (CCTK_Equals(lapseinit,"schwarz")==1) then 
         print *,"Initialise with schwarzschild lapse"
         alp = (2.*r - mass)/(2.*r+mass)
      endif
      
c     time symmetric initial slice
      hxx = zero 
      hxy = zero
      hxz = zero
      hyy = zero
      hyz = zero
      hzz = zero

      return
      end subroutine Schwarzschild
