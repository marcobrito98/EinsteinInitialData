/*@@
  @file      brilldata.F
  @date
  @author    Carsten Gundlach, Miguel Alcubierre.
  @desc
             Construct Brill wave q function.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"

      function brillq(rho,z,phi)

c     Calculates the function q that appear in the conformal
c     metric for Brill waves:
c
c     ds^2  =  psi^4 ( e^(2q) (drho^2 + dz^2) + rho^2 dphi^2 )
c
c     There are three different choices for the form of q depending
c     on the value of the parameter "brill_q":
c
c     brill_q = 0:
c                                           2
c                           2+b  - (rho-rho0)      2     2
c                  q = a rho    e            (z/sz)  /  r
c
c
c     brill_q = 1:  (includes Eppleys form as special case)
c
c
c                                  b               2    2      2  c/2
c                  q = a (rho/srho)  /  { 1 + [ ( r - r0 ) / sr  ]    }
c
c
c     brill_q = 2:  (includes Holz et al form as special case)
c
c                                            2    2      2  c/2 
c                                  b  - [ ( r - r0 ) / sr  ]
c                  q = a (rho/srho)  e
c
c
c     If we want 3D initial data, the function q is multiplied by an
c     additional factor:
c
c                           m    2                              m
c         q -> q [ 1 + d rho  cos  (n (phi)+phi0 ) / ( 1 + e rho ) ]

      implicit none

      DECLARE_CCTK_PARAMETERS

      logical firstcall

      integer qtype
      integer b,c

      CCTK_REAL brillq,rho,z,phi
      CCTK_REAL a,r0,sr,rho0,srho,sz
      CCTK_REAL d,e,m,n,phi0

      data firstcall /.true./

      save firstcall,qtype
      save a,b,c,r0,sr,rho0,srho,sz
      save d,e,m,n,phi0

c     Get parameters at first call.

      if (firstcall) then

         qtype = brill_q

         a = brill_a
         b = brill_b
         c = brill_c

         r0   = brill_r0
         sr   = brill_sr
         rho0 = brill_rho0
         srho = brill_srho
         sz   = brill_sz

         if (axisym.eq.0) then
            d = brill_d
            e = brill_e
            m = brill_m
            n = brill_n
            phi0 = brill_phi0
         end if

         firstcall = .false.

      end if

c     Calculate q(rho,z) from a choice of forms.

c     brill_q = 0.

      if (qtype.eq.0) then

         brillq = a*dabs(rho)**(2 + b)
     .          / dexp((rho - rho0)**2)/(rho**2 + z**2)

         if (sz.ne.0.0D0) then
            brillq = brillq/dexp((z/sz)**2)
         end if

      else if (qtype.eq.1) then

c     brill_q = 1.  This includes Eppleys choice of q. 
c     But note that q(Eppley) = 2q(Cactus).

         brillq = a*(rho/srho)**b 
     .          / ( 1.0D0 + ((rho**2 + z**2 - r0**2)/sr**2)**(c/2) )

      else if (qtype.eq.2) then

c     brill_q = 2.   This includes my (Carstens) notion of what a
c     smooth "pure quadrupole" q should be, plus the choice of
c     Holz et al.

         brillq = a*(rho/srho)**b
     .          / dexp(((rho**2 + z**2 - r0**2)/sr**2)**(c/2))

         print *,a,b,c,r0,sr,srho,z,rho,brillq

      else

c     Unknown type for q function.

         call CCTK_WARN(0,"Unknown type of Brill function q")

      end if

c     If desired, multiply with a phi dependent factor.

      if (axisym.eq.0) then
         brillq = brillq*(1.0D0 + d*rho**m*cos(n*(phi-phi0))**2
     .          / (1.0D0 + e*rho**m))
      end if

      return
      end



