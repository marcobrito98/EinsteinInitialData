#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      function brillq(rho,z,phi)

C     Authors:  Carsten Gundlach, Miguel Alcubierre.
C
C     Calculates the function q that appear in the conformal
C     metric for Brill waves:
C
C     ds^2  =  psi^4 ( e^(2q) (drho^2 + dz^2) + rho^2 dphi^2 )
C
C     There are three different choices for the form of q depending
C     on the value of the parameter "brill_q":
C
C     brill_q = 0:
C                                           2
C                           2+b  - (rho-rho0)      2     2
C                  q = a rho    e            (z/sz)  /  r
C
C
C     brill_q = 1:  (includes Eppleys form as special case)
C
C
C                                  b               2    2      2  c/2
C                  q = a (rho/srho)  /  { 1 + [ ( r - r0 ) / sr  ]    }
C
C
C     brill_q = 2:  (includes Holz et al form as special case)
C
C                                            2    2      2  c/2 
C                                  b  - [ ( r - r0 ) / sr  ]
C                  q = a (rho/srho)  e
C
C
C     If we want 3D initial data, the function q is multiplied by an
C     additional factor:
C
C                           m    2                              m
C         q -> q [ 1 + d rho  cos  (n (phi)+phi0 ) / ( 1 + e rho ) ]


      implicit none

      DECLARE_CCTK_PARAMETERS

      logical firstcall

      integer CCTK_Equals
      integer qtype

      real*8 brillq,rho,z,phi
      real*8 a,b,c,r0,sr,rho0,srho,brill_sz
      real*8 d,e,m,n,phi0

      data firstcall /.true./

      save firstcall,qtype,a,b,c,r0,sr,rho0,srho,brill_sz
      save d,e,m,n,phi0

C     Get parameters at first call.

      if (firstcall) then

         qtype = brill_q

         a = brill_a
         b = brill_b
         c = brill_c

         r0   = brill_r0
         sr   = brill_sr
         rho0 = brill_rho0
         srho = brill_srho
         brill_sz = brill_sz

         if (axisym.eq.0) then
            d = brill_d
            e = brill_e
            m = brill_m
            n = brill_n
            phi0 = brill_phi0
         end if

         firstcall = .false.

      end if

      if (rho.lt.0) call CCTK_WARN(1,"Warning: negative rho in brillq:")

C     Calculate q(rho,z) from a choice of forms. Type 0, 1 and 2 are those
C     of Bruegmann.

C     brill_q = 0.

      if (qtype .eq. 0) then

         brillq = a * rho**(2.d0+b)
     $        / dexp((rho - rho0)**2) / (rho**2 + z**2) 
         if (brill_sz .ne. 0.d0) then
            brillq = brillq / exp((z/brill_sz)**2)
         end if

      else if (qtype .eq. 1) then

C     brill_q = 1.  This includes Eppleys choice of q. 
C     But note that q(Eppley) = 2q(Cactus).

         brillq = a * (rho/srho)**b 
     $        / ( 1.d0 + ((rho**2 + z**2 - r0**2) / sr**2)**(c/2) )

      else if (qtype .eq. 2) then

C     brill_q = 2.   This includes my (Carstens) notion of what a
C     smooth "pure quadrupole" q should be.

         brillq = a * (rho/srho)**b 
     $        / dexp( ((rho**2 + z**2 - r0**2) / sr**2)**(c/2) )

      else

C     Unknown type for q function.

         call CCTK_WARN(0,"Unknown type of Brill function q.")

      end if

C     If desired, multiply with phi dependent factor.

      if (axisym.eq.0) then
         brillq = brillq*(1.0D0 + d*rho**m*cos(n*(phi-phi0))**2
     .          /(1.0D0 + e*rho**m))
      end if

      return
      end



