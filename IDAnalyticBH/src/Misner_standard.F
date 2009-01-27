 /*@@
   @file      Misner_standard.F
   @date      March 1997
   @author    Joan Masso
   @desc 
      Set up initial data for two Misner black holes
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_arguments.h"
#include "cctk_parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

 /*@@
   @routine    Misner_standard
   @date       
   @author     Joan Masso, Ed Seidel
   @desc 
            Initialize the metric with a time symmetrical 
            black hole spacetime containing
            two axially symmetric misner black holes with a 
            mass/length parameter mu. The mass is computed.
            The spacetime line element has the form:
                 $$ ds^2 = -dt^2 + \Psi^4 (dx^2+dy^2+dz^2) $$
            and only $\Psi$ differs.
            (Conformal factor from Karen Camarda)
   @enddesc 
   @calls     
   @history 
 
   @endhistory 

   @par      mu
   @pdesc    Misner parameter. 
   @ptype    real
   @pcomment Values less than 1.8 do not really correspond to two 
            black holes, as there is an initial single event horizon 
            surrounding the throats. So, with low values of mu we also
            have distorted single black holes.
   @endpar
 
   @par      nmax
   @pdesc    Summation limit for the misner series in the 'twobh' case.
   @ptype    integer
   @endpar

@@*/
      subroutine Misner_standard(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer :: n
      INTEGER :: CCTK_Equals
      CCTK_REAL :: zero, one, three

      zero  = 0.d0
      one   = 1.d0
      three = 3.d0

c     Initialize so we can accumulate
c     -------------------------------
      psi = one
      if (use_conformal_derivs == 1) then
         psix = zero
         psiy = zero
         psiz = zero
         psixx = zero
         psixy = zero
         psixz = zero
         psiyy = zero
         psiyz = zero
         psizz = zero
      end if

      do n = nmax,1,-1
         csch = one/sinh(mu*n)
         coth = one/tanh(mu*n)
         r1 = sqrt(x**2+y**2+(z+coth)**2)
         r2 = sqrt(x**2+y**2+(z-coth)**2)
         psi = psi + csch*(one/r1 + one/r2)

         if (use_conformal_derivs == 1) then
            psix = psix + (-(x/r2**3)-x/r1**3)*csch
            psiy = psiy + (-(y/r2**3)-y/r1**3)*csch
            psiz = psiz + (-((z-coth)/r2**3)-(z + coth)/r1**3)*csch
            psixx = psixx +
     &           ((three*x**2)/r2**5-one/r2**3+(three*x**2)/r1**5-one/r1**3)*csch
            psixy = psixy + ((three*x*y)/r2**5+(three*x*y)/r1**5)*csch
            psixz = psixz + ((three*x*(z-coth))/r2**5
     &           +(three*x*(z+coth))/r1**5)*csch
            psiyy = psiyy +
     &           ((three*y**2)/r2**5-one/r2**3
     &           +(three*y**2)/r1**5-one/r1**3)*csch
            psiyz = psiyz + ((three*y*(z-coth))/r2**5
     &           +(three*y*(z+coth))/r1**5)*csch
            psizz = psizz + (-one/r2**3+(three*(z-coth)**2)/r2**5+
     &           (three*(z+coth)**2)/r1**5-one/r1**3)*csch

         end if

      end do

c     Cactus convention
c     -----------------
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

c     compute the ADM mass
c     -------------------- 
      mass = zero
      do n = 1,nmax
         mass = mass + 4./sinh(n*mu)
      enddo
      write(*,*) "   >>ADM mass",mass

c     Should initialize lapse to Cadez value if possible
c     --------------------------------------------------
      if (CCTK_Equals(initial_lapse,"cadez") == 1) then 
         write (*,*)"Initialise with cadez lapse"
         alp = one
         do n = 1,nmax
            coth = one/tanh(mu*n)
            r1 = sqrt(x**2+y**2+(z+coth)**2)
            r2 = sqrt(x**2+y**2+(z-coth)**2)
            alp = alp + (-1.)**n * one/sinh(mu*n)*(one/r1 + one/r2)
         enddo
         alp = alp/psi
      endif
         
      
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
