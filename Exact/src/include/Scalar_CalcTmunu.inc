C Here we added the matter variables for several of the metrics
C you can find in "src" directory through the components of the
C stress-eergy tensor. Being different in specific cases it is
C necessary to check first which metric is running.  
C Author : Dumitru Vulcanov (Timisoara, Romania)
C $Header$

C Varianta cu un singur param. (raza initiala) pt Rob-Walker

#include "param_defs.inc"

c
c FIXME:
c If we could be certain that this code were always compiled as Fortran 90,
c the decode here could be done with a case statement, which ought to give
c better performance than the if-else chain we use now.  But in practice
c Cactus has enough other overheads that this is not bothering with...
c

c
c ***** KLUDGE *****
c
c This code is #include-d into various evolution thorns, and alas does not
c have direct access to thorn Exact parameters.  Instead, this code must
c use the restricted-grid-scalar copies of the parameters.  In practice,
c this means changing "__" to "___" in all parameter names (but not in the
c #define constants in "param_defs.inc").  See the comments in param.ccl
c for further information on this.
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Minkowski spacetime cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if     (decoded_exact_model .eq. EXACT__Minkowski) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Minkowski_shift) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Minkowski_funny) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Minkowski_gauge_wave) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Minkowski_conf_wave) then
c        no stress-energy tensor in this model

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc black hole spacetimes cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       elseif (decoded_exact_model .eq. EXACT__Schwarzschild_EF) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Schwarzschild_PG) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Schwarzschild_Novikov) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Kerr_BoyerLindquist) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Kerr_KerrSchild) then
c        no stress-energy tensor in this model

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c Schwarzschild-Lemaitre spacetime
c (Schwarzschild black hole with cosmological constant)
c
       elseif (decoded_exact_model .eq. EXACT__Schwarzschild_Lemaitre) then
         razsch2=x(i,j,k)*x(i,j,k)+y(i,j,k)*y(i,j,k)+z(i,j,k)*z(i,j,k)
         coefsch=-Schwarzschild_Lemaitre___Lambda/(8.0D0*pi)
         pppsch=1.0D0-2.0D0*Schwarzschild_Lemaitre___mass/sqrt(razsch2)
     &               -Schwarzschild_Lemaitre___Lambda*razsch2/3.0D0
         unusch=(1.0D0-pppsch)/pppsch/razsch2

         Ttt =  Ttt-coefsch*pppsch
         Ttx =  Ttx
         Tty =  Tty
         Ttz =  Ttz
         Txx  = Txx+coefsch*(1.0D0+x(i,j,k)*x(i,j,k)*unusch)
         Tyy  = Tyy+coefsch*(1.0D0+y(i,j,k)*y(i,j,k)*unusch)
         Tzz  = Tzz+coefsch*(1.0D0+z(i,j,k)*z(i,j,k)*unusch)
         Txy  = Txy+coefsch*x(i,j,k)*y(i,j,k)*unusch
         Txz  = Txz+coefsch*x(i,j,k)*z(i,j,k)*unusch
         Tyz  = Tyz+coefsch*y(i,j,k)*z(i,j,k)*unusch

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       elseif (decoded_exact_model .eq. EXACT__multi_BH) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Alvi) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Thorne_fakebinary) then
c        no stress-energy tensor in this model

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc cosmological spacetimes cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c Lemaitre spacetime
c
        elseif (decoded_exact_model .eq. EXACT__Lemaitre) then
          unu1 = sqrt(3.0D0*Lemaitre___Lambda)
     &           * CCTK_time * (Lemaitre___kappa+1.0D0) / (2.0D0)
          raz  = Lemaitre___R0*(cosh(unu1)
     &                         + sqrt(1.0D0+8.0D0*pi*Lemaitre___epsilon0
     &                                               /Lemaitre___Lambda)
     &                           *sinh(unu1))
     &                        **(2.0D0/(3.0D0*Lemaitre___kappa+3.0D0)) 
          raz2 = raz*raz
          treiori = -Lemaitre___Lambda*raz2/8.0D0/pi
     &           +Lemaitre___epsilon0*Lemaitre___kappa
     &                              *raz**(-3.0D0*Lemaitre___kappa-1.0D0)

          Ttt = Ttt + Lemaitre___Lambda/8.0D0/pi
     &              + Lemaitre___epsilon0*raz**(-3.0D0*(Lemaitre___kappa+1.0D0))
          Txx = Txx + treiori
          Tyy = Tyy + treiori
          Tzz = Tzz + treiori
          Txy = Txy
          Tyz = Tyz
          Txz = Txz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc this metric doesnt work and has been moved to ../../archive/
ccc
ccc Robertson-Walker spacetime
ccc
cc        elseif (decoded_exact_model .eq. EXACT__Robertson_Walker) then
cc          rr2  = x(i,j,k)*x(i,j,k)+y(i,j,k)*y(i,j,k)+z(i,j,k)*z(i,j,k)
cc
cc          if (Robertson_Walker___pressure .gt. 0) then
cc            aha1  = Robertson_Walker___k * (Robertson_Walker___R0**2)
cc     &                                  / (8.0D0*pi*(raza(i,j,k)**2))
cc            aha2 = Robertson_Walker___k/(1.0D0 - Robertson_Walker___k*rr2) 
cc
cc            Ttt =  Ttt + 3.0D0*aha1/(raza(i,j,k)*raza(i,j,k))
cc            Txx =  Txx + aha1*(1.0D0 + aha2*x(i,j,k)*x(i,j,k))
cc            Tyy =  Tyy + aha1*(1.0D0 + aha2*y(i,j,k)*y(i,j,k))
cc            Tzz =  Tzz + aha1*(1.0D0 + aha2*z(i,j,k)*z(i,j,k))
cc            Txy =  Txy + aha1*aha2*x(i,j,k)*y(i,j,k)
cc            Txz =  Txz + aha1*aha2*x(i,j,k)*z(i,j,k)
cc            Tyz =  Tyz + aha1*aha2*y(i,j,k)*y(i,j,k)
cc          else
cc            Ttt = Ttt+Robertson_Walker___rho * (Robertson_Walker___R0**3)
cc     &                                      / (raza(i,j,k)**3)
cc            Txx = Txx
cc            Tyy = Tyy
cc            Tzz = Tzz
cc            Txy = Txy
cc            Txz = Txz
cc            Tyz = Tyz
cc          endif
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c de Sitter spacetime
c
        elseif (decoded_exact_model .eq. EXACT__de_Sitter) then
          Ttt =  Ttt + 1.0D0/6.0D0/pi/(CCTK_time**2)
          Ttx  = Ttx
          Tty  = Tty
          Ttz  = Ttz
          Txx  = Txx 
          Tyy  = Tyy 
          Tzz  = Tzz 
          Txy  = Txy
          Txz  = Txz
          Tyz  = Tyz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c de Sitter spacetime with cosmological constant
c
        elseif (decoded_exact_model .eq. EXACT__de_Sitter_Lambda) then
          aaaa = de_Sitter_Lambda___scale/(8.0D0*pi)
          bbbb = aaaa*exp(2.0D0*sqrt(de_Sitter_Lambda___scale/3.0D0)*CCTK_time)

          Ttt  = Ttt + aaaa
          Ttx  = Ttx 
          Tty  = Tty
          Ttz  = Ttz
          Txx  = Txx - bbbb
          Tyy  = Tyy - bbbb
          Tzz  = Tzz - bbbb
          Txy  = Txy
          Txz  = Txz
          Tyz  = Tyz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c anti-de Sitter spacetime with cosmological constant
c

        elseif (decoded_exact_model .eq. EXACT__anti_de_Sitter_Lambda) then
          aaaa1 = anti_de_Sitter_Lambda___scale/(8.0D0*pi)
          bbbb1 = aaaa1*exp(2.0D0*sqrt(-anti_de_Sitter_Lambda___scale/3.0D0)
     &                           *x(i,j,k))
        
          Ttt =  Ttt + bbbb1
          Ttx  = Ttx 
          Tty  = Tty
          Ttz  = Ttz
          Txx  = Txx - aaaa1
          Tyy  = Tyy - bbbb1
          Tzz  = Tzz - bbbb1
          Txy  = Txy
          Txz  = Txz
          Tyz  = Tyz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       elseif (decoded_exact_model .eq. EXACT__Bianchi_I) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__Goedel) then
c        no stress-energy tensor in this model

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c Bertotti spacetime
c
        elseif (decoded_exact_model .eq. EXACT__Bertotti) then
          bass  = Bertotti___Lambda/(8.0D0*pi)
          term1 = bass*exp(2.0D0*sqrt(-Bertotti___Lambda)*x(i,j,k))
          term2 = bass*exp(2.0D0*sqrt(-Bertotti___Lambda)*z(i,j,k))

          Ttt =  Ttt + term1
          Ttx  = Ttx 
          Tty  = Tty
          Ttz  = Ttz
          Txx  = Txx - bass
          Tyy  = Tyy - term2
          Tzz  = Tzz - bass
          Txy  = Txy
          Txz  = Txz
          Tyz  = Tyz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c Kasner-like spacetime
c
        elseif (decoded_exact_model .eq. EXACT__Kasner_like) then

          kkkk=Kasner_like___q*(2.0D0-3.0D0*Kasner_like___q)
     &                        /(8.0D0*pi*(CCTK_time**2))

          Ttt  = Ttt + kkkk
          Ttx  = Ttx
          Tty  = Tty
          Ttz  = Ttz
          Txx  = Txx + kkkk*CCTK_time**(2.0D0*Kasner_like___q)
          Tyy  = Tyy + kkkk*CCTK_time**(2.0D0*Kasner_like___q)
          Tzz  = Tzz + kkkk*CCTK_time**(2.0D0-4.0D0*Kasner_like___q)
          Txy  = Txy
          Txz  = Txz
          Tyz  = Tyz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       elseif (decoded_exact_model .eq. EXACT__Kasner_axisymmetric) then
c        no stress-energy tensor in this model

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c generalized Kasner spacetime
c
        elseif (decoded_exact_model .eq. EXACT__Kasner_generalized) then

        kkkk1 = (  Kasner_generalized___p1 - Kasner_generalized___p1**2
     &           + Kasner_generalized___p2 - Kasner_generalized___p2**2
     &           - Kasner_generalized___p1*Kasner_generalized___p2 )
     &          / (8.0D0*pi*(CCTK_time**2))

        Ttt  = Ttt + kkkk1
        Ttx  = Ttx
        Tty  = Tty
        Ttz  = Ttz
        Txx  = Txx+kkkk1*CCTK_time**(2.0D0*Kasner_generalized___p1)
        Tyy  = Tyy+kkkk1*CCTK_time**(2.0D0*Kasner_generalized___p2)
        Tzz  = Tzz+kkkk1*CCTK_time**(2.0D0-2.0D0*Kasner_generalized___p1
     &                                    -2.0D0*Kasner_generalized___p2)
        Txy  = Txy
        Txz  = Txz
        Tyz  = Tyz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (decoded_exact_model .eq. EXACT__Gowdy_wave) then
c        no stress-energy tensor in this model


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       elseif (decoded_exact_model .eq. EXACT__Milne) then
c        no stress-energy tensor in this model

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc miscellaneous spacetimes ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       elseif (decoded_exact_model .eq. EXACT__boost_rotation_symmetric) then
c        no stress-energy tensor in this model
       elseif (decoded_exact_model .eq. EXACT__bowl) then
c        no stress-energy tensor in this model

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c constant density star
c
        elseif (decoded_exact_model .eq. EXACT__constant_density_star) then
          rr  = sqrt(x(i,j,k)*x(i,j,k)+y(i,j,k)*y(i,j,k)+
     &                z(i,j,k)*z(i,j,k))
          star_m = constant_density_star___mass
          star_r = constant_density_star___radius

          r3 = star_r**3
          if (rr.le.star_r) then
            unu = 3.0D0*sqrt(1.0D0-2.0D0*star_m/star_r)
            doi = sqrt(1.0D0-2.0D0*star_m*rr*rr/r3)
            trei= star_m*(unu-3.0D0*doi)/(2*pi*(unu-doi)*r3)
            Ttt = Ttt + 3.0D0*star_m*
     &        (5.0D0-9.0D0*star_m/star_r - unu*doi
     &        -star_m*rr*rr/r3)/(8.0D0*pi*r3)
            Txx = Txx -trei*(1.0D0+2.0D0*star_m*x(i,j,k)*x(i,j,k)/
     &      (doi*doi*r3))/2.0D0
            Tyy = Tyy -trei*(1.0D0+2.0D0*star_m*y(i,j,k)*y(i,j,k)/
     &      (doi*doi*r3))/2.0D0
            Tzz = Tzz -trei*(1.0D0+2.0D0*star_m*z(i,j,k)*z(i,j,k)/
     &      (doi*doi*r3))/2.0D0
            Txy = Txy -trei*star_m*x(i,j,k)*y(i,j,k)/(doi*doi*r3)
            Tyz = Tyz -trei*star_m*y(i,j,k)*z(i,j,k)/(doi*doi*r3)
            Txz = Txz -trei*star_m*x(i,j,k)*z(i,j,k)/(doi*doi*r3)
          endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        else
          call CCTK_WARN(0,"Unknown value of Exact::decoded_exact_model")
        endif 
