#include "cctk.h"
#include "declare_arguments.h"
#include "declare_parameters.h"

/*@@
  @file     teukwaves.F
  @date     Jan 96
  @author   Joan Masso + Malcolm Tobias
  @desc
            Routines for the teukolsky waves initialization.
  @enddesc

@@*/
        
/*@@
  @routine  teukwaves
  @date     Jan 96
  @author   Malcolm Tobias + Joan Masso 
  @desc 
            Generates the analytic solutions of the
            metric components for linearized quadrupolar (teukolsky) waves.
            More changes for cactus port: all 3d arrays are converted
            to point calculations. What a pain!!!

            Three and a half years later....
            Converted to Cactus4.0

  @enddesc
@@*/

      subroutine teukwaves(CCTK_FARGUMENTS)
      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_PARAMETERS

      
      REAL amp,m,ra,pi
      REAL wave,wave2,wave3,wave4,wave5,wave6,wave7,wave8
      REAL wavepulse
      INTEGER ipacket,iparity
      INTEGER ingoing,outgoing
      REAL kappa
c     point values of x.y z,r
      REAL xp,yp,zp,rp

c     spherical metric
      REAL grr,grt,grp,gtt,gtp,gpp
      REAL hrr,hrt,hrp,htt,htp,hpp
      
c     special coefficeints needed for teukolsky waves
c     All these were 3d arrays in old versions of Newage,G,etc..
c     uffff....
      REAL teuk_tp,teuk_tn
      REAL tp2,tn2,fp,fn,fpa,fna,fpb,fnb,fpc,fnc,fpd,fnd,fpe,fne
      REAL ca,cb,cc,ck,cl,frr,frt,frp,ftt1,ftt2,ftp,fpp1,fpp2
      REAL cadot,cbdot,ccdot,ckdot,cldot
      REAL drt,drp,dtt,dtp,dpp,sint,cost,sinp,cosp,tmp

c     from old spheretocart
      REAL drx,dry,drz,dtx,dty,dtz,dpx,dpy,dpz
      
      INTEGER i,j,k

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      pi    = 3.14159265358979
      kappa = sqrt(143.0d0/pi)/12288.0d0

      if (wavecenter.ne.0.) then
         CCTK_Warn(0,"Teukwaves: wavecenter must be zero for time symmetry")
      endif
 
c     CACTUS 4.0: all REAL/INTEGER parameter are available as variables and 
c     need not be derived from the database. Nevertheless I will assign 
c     the parameter variables to the variables of CAC3.0, to allow for 
c     easy portation (wavepulse stays wavepulse btw). (GERD)
      
      am  = amplitude
      m   = mvalue
      wave= wavelength

c     to be consistent with Evans & Abrahams, rescale the amplitude
c     by kappa, because they use that in the packet and
c     4/3 because they use I not Teuk. F  
      if (CCTK_Equals(packet,'eppley')) then
         ipacket = 1
         write(*,*)'Teuk.Packet= eppley' 
      elseif (CCTK_Equals(packet,'evans')) then
         ipacket = 2
         write(*,*)'Teuk.Packet= evans'
c     to be consistent with Evans & Abrahams, rescale the amplitude
c     by kappa, because they use that in the packet and
c     4/3 because they use I not Teuk. F
         amp = amp*(4.0d0/3.0d0)*kappa
      elseif (CCTK_Equals(packet,'square')) then
         ipacket = 3
         write(*,*)'Teuk.Packet= square'
      else
         write(*,*)'Teukwaves: unsupported packet'
         STOP
      endif      

      
      write(*,*)'amplitude  = ',amp
      write(*,*)'     m     = ',m
      write(*,*)'wavecenter = ',ra
      write(*,*)'wavelength = ',wave
      
      if (CCTK_Equals(parity,'even')) then    
         iparity = 0
      else (CONTAINS(parity,'odd')) then 
         iparity = 1
      endif

       if(CCTK_Equals(wavegoing,'in')) then
          ingoing = 1
       endif

      if(CCTK_Equals(wavegoing,'out')) then
         outgoing = 1
      endif

      if(CCTK_Equals(wavegoing,'both')) then
         ingoing  = 1
         outgoing = 1
      endif

      if (ipacket.eq.1) then
         if (ingoing.eq.0.or.outgoing.eq.0) then
            write(*,*)'WARNING: Epply packet is only non-singular at'
            write(*,*)'the origin for ingoing-outgoing combination'
            write(*,*)'of waves.  Expect problems ;-) '
         endif
      endif
      
      wave2 = wave*wave
      wave3 = wave2*wave
      wave4 = wave3*wave
      wave5 = wave4*wave
      wave6 = wave5*wave
      wave7 = wave6*wave
      wave8 = wave7*wave
         

c  ****************************************   
c     initial data for teukolsky waves
c  ****************************************
      do k=1,sh(3)
         do j=1,sh(2)
            do i=1,sh(1)
               
               xp = x(i,j,k)
               yp = y(i,j,k)
               zp = z(i,j,k)
               rp = r(i,j,k)
               
               teuk_tp = (time+rp-ra)
               teuk_tn = (time-rp+ra)
               
               tp2 = teuk_tp**2
               tn2 = teuk_tn**2

c initialize everything to zero
               fp = 0.0d0
               fpa = 0.0d0
               fpb = 0.0d0
               fpc = 0.0d0
               fpd = 0.0d0
               fpe = 0.0d0

               fn = 0.0d0
               fna = 0.0d0
               fnb = 0.0d0
               fnc = 0.0d0
               fnd = 0.0d0
               fne = 0.0d0
               
c     eppley package --> x exp(-x^2)  
               if (ipacket.eq.1) then
c     to keep the analytic solution valid for the eppley packet:   
                  if (ra.ne.0.and.ingoing.eq.1.and.outgoing.eq.1) then
                     teuk_tn = (time-rp-ra)
                     tn2 = teuk_tn**2
                  endif

                  if (ingoing.eq.1) then

                     tmp = exp(-tp2)
                     fp = amp * teuk_tp*tmp
                     fpa = amp * (1 - 2*tp2)*tmp
                     fpb = amp * teuk_tp*(4*tp2 - 6)*tmp
                     fpc = amp * (24*tp2 - 8*tp2*tp2 - 6)*tmp
                     fpd = amp * teuk_tp*(60 - 80*tp2 + 16*tp2*tp2)*tmp
                     fpe = amp * (
     &                    (60 - 360*tp2 + 240*tp2*tp2 - 32*tp2*tp2*tp2)*tmp )

                     if (outgoing.eq.0.and.(teuk_tp+wavepulse).le.0) 
     &                    then

                        fp  = 0.0d0
                        fpa = 0.0d0
                        fpb = 0.0d0
                        fpc = 0.0d0
                        fpd = 0.0d0
                        fpe = 0.0d0

                     endif

                  endif

                  if (outgoing.eq.1) then

                     tmp = exp(-tn2)
                     fn = amp * teuk_tn*tmp
                     fna = amp * (1 - 2*tn2)*tmp
                     fnb = amp * teuk_tn*(4*tn2 - 6)*tmp
                     fnc = amp * (24*tn2 - 8*tn2*tn2 -6)*tmp
                     fnd = amp * teuk_tn*(60 - 80*tn2 + 16*tn2*tn2)*tmp
                     fne = amp * (
     &                    (60 - 360*tn2 + 240*tn2*tn2 - 32*tn2*tn2*tn2)*tmp )

                  endif
               endif

      
c     evans package  --> w^4(1-x^2/w^2)^6
               if(ipacket.eq.2) then

                  if (ingoing.eq.1) then

                     tmp = 1.0d0 - tp2/wave2

                     if (abs(teuk_tp).lt.wave) then

                        fp  =  amp*wave5*tmp**6
                        fpa =  -12.0d0*amp*wave3*teuk_tp*tmp**5
                        fpb =  -12.0d0*amp*wave3*(
     &                       tmp**5 - 10.0d0*tp2/wave2*tmp**4 )
                        fpc =  120.0d0*amp*wave3*(
     &                        3.0d0*teuk_tp/wave2*tmp**4 
     &                       -8.0d0*teuk_tp*tp2/wave4*tmp**3 )
                        fpd = 360.0d0*amp*wave3*(
     &                        1.0d0/wave2*tmp**4 
     &                       -16.0d0*tp2/wave4*tmp**3
     &                       +16.0d0*tp2*tp2/wave6*tmp**2 )
                        fpe = 2880.0d0*amp*wave3*(
     &                       -5.0d0*teuk_tp/wave4*tmp**3
     &                       +20.0d0*teuk_tp*tp2/wave6*tmp**2
     &                       -8.0d0*teuk_tp*tp2*tp2/wave8*tmp )

                     else
                        fp=0.
                        fpa=0.
                        fpb=0.
                        fpc=0.
                        fpd=0.
                        fpe=0.
                     endif

                  endif

                  if (outgoing.eq.1) then

                     tmp = 1.0d0 - tn2/wave2
         
                     if (abs(teuk_tn).lt.wave) then

                        fn  =  amp*wave5*tmp**6
                        fna =  -12.0d0*amp*wave3*teuk_tn*tmp**5
                        fnb =  -12.0d0*amp*wave3*(
     &                       tmp**5 - 10.0d0*tn2/wave2*tmp**4 )
                        fnc =  120.0d0*amp*wave3*(
     &                        3.0d0*teuk_tn/wave2*tmp**4 
     &                       -8.0d0*teuk_tn*tn2/wave4*tmp**3 )
                        fnd = 360.0d0*amp*wave3*(
     &                        1.0d0/wave2*tmp**4 
     &                       -16.0d0*tn2/wave4*tmp**3
     &                       +16.0d0*tn2*tn2/wave6*tmp**2 )
                        fne = 2880.0d0*amp*wave3*(
     &                       -5.0d0*teuk_tn/wave4*tmp**3
     &                       +20.0d0*teuk_tn*tn2/wave6*tmp**2
     &                       -8.0d0*teuk_tn*tn2*tn2/wave8*tmp )

                     else
                        fn=0.
                        fna=0.
                        fnb=0.
                        fnc=0.
                        fnd=0.
                        fne=0.
                     endif

                  endif

               endif
      
c     square package --> (1-x^2/w^2)^2
               if(ipacket.eq.3) then 

                  write(*,*)'WARNING: need to calculate fpe,fne'
                  write(*,*)'WARNING: need conditionals for in out waves'

                  if (abs(teuk_tp).lt.wave) then
                     fp = amp * (1-(tp2/wave2))**2
                     fpa = amp * (-4*teuk_tp)*(1-tp2/wave2)/wave2
                     fpb = amp * (8*tp2/wave4-4*(1-tp2/wave2)/wave2)
                     fpc = amp * 24*teuk_tp/wave4
                     fpd = amp * 24/wave4
                  else
                     fp=0.
                     fpa=0.
                     fpb=0.
                     fpc=0.
                     fpd=0.
                  endif
         
                  if (abs(teuk_tn).lt.wave) then
                     fn = amp * (1-(tn2/wave2))**2
                     fna = amp * (-4*teuk_tn)*(1-tn2/wave2)/wave2
                     fnb = amp * (8*tn2/wave4-4*(1-tn2/wave2)/wave2)
                     fnc = amp * 24*teuk_tn/wave4
                     fnd = amp * 24/wave4  
                  else
                     fn=0.
                     fna=0.
                     fnb=0.
                     fnc=0.
                     fnd=0.
                  endif
               endif
      
c     
c     coefficients
c     
               ca = 3*( (fpb-fnb)/rp**3 -3*(fna+fpa)/rp**4 +3*(fp-fn)/rp**5 )
               cb = -( -(fnc+fpc)/rp**2 +3*(fpb-fnb)/rp**3 -6*(fna+fpa)/rp**4
     &              +6*(fp-fn)/rp**5 ) 
               cc = ( (fpd-fnd)/rp -2*(fnc+fpc)/rp**2 +9*(fpb-fnb)/rp**3
     &              -21*(fna+fpa)/rp**4 +21*(fp-fn)/rp**5 )/4
               ck = (fpb-fnb)/rp**2-3*(fpa+fna)/rp**3+3*(fp-fn)/rp**4
               cl = -(fpc+fnc)/rp+2*(fpb-fnb)/rp**2-3*(fpa+fna)/rp**3+
     &              3*(fp-fn)/rp**4
               
c     If we only have an outgoing wave, we should flip the sign
c     of the above coefficients.  That is because the above coefficients
c     assume an ingoing and outgoing wave packet of the form:
c     wave = ingoing - outgoing.  For a single wave, we want ingoing
c     and outgoing to have the same sign.

               if(ingoing.eq.0.and.outgoing.eq.1) then
                  ca=-ca
                  cb=-cb
                  cc=-cc
                  ck=-ck
                  cl=-cl
               endif
               
               sint = (xp**2+yp**2)**0.5/rp
               cost = zp/rp
               sinp = yp/(xp**2+yp**2)**0.5
               cosp = xp/(xp**2+yp**2)**0.5
               
c     mvalue
               if (m.eq.0) then
                  frr = 2-3*sint**2
                  frt = -3*sint*cost
                  frp = 0.
                  ftt1 = 3*sint**2
                  ftt2 = -1.
                  ftp = 0.
                  fpp1 = -ftt1
                  fpp2 = 3*sint**2-1
                  drt = 0.
                  drp = -4*cost*sint
                  dtt = 0.
                  dtp = -sint**2
                  dpp = 0.
               elseif (m.eq.-1) then
                  frr = 2*sint*cost*sinp
                  frt = (cost**2-sint**2)*sinp 
                  frp = cost*cosp 
                  ftt1 = -frr
                  ftt2 = 0.
                  ftp = sint*cosp 
                  fpp1 = -ftt1
                  fpp2 = ftt1
                  drt = -2*cost*sinp
                  drp = -2*(cost**2-sint**2)*cosp
                  dtt = -sint*sinp
                  dtp = -cost*sint*cosp
                  dpp = sint*sinp
               elseif (m.eq.-2) then
                  frr = sint**2*2*sinp*cosp
                  frt = sint*cost*2*sinp*cosp
                  frp = sint*(1-2*sinp**2)
                  ftt1 = (1+cost**2)*2*sinp*cosp
                  ftt2 = -2*sinp*cosp
                  fpp1 = -ftt1
                  fpp2 = cost**2*2*sinp*cosp
                  drt = 8*sint*sinp*cosp
                  drp = 4*sint*cost*(1-2*sinp**2)
                  dtt = -4*cost*sinp*cosp
                  dtp = (2-sint**2)*(2*sinp**2-1)
                  dpp = 2*cost*2*sinp*cosp
               elseif (m.eq.1) then
                  frr = 2*sint*cost*cosp
                  frt = (cost**2-sint**2)*cosp
                  frp = -cost*sinp
                  ftt1 = -2*sint*cost*cosp
                  ftt2 = 0.
                  ftp = -sint*sinp
                  fpp1 = -ftt1
                  fpp2 = ftt1
                  drt = -2*cost*cosp
                  drp = 2*(cost**2-sint**2)*sinp
                  dtt = -sint*cosp
                  dtp = cost*sint*sinp
                  dpp = sint*cosp
               elseif (m.eq.2) then
                  frr = sint**2*(1-2*sinp**2)
                  frt = sint*cost*(1-2*sinp**2)
                  frp = -sint*2*sinp*cosp
                  ftt1 = (1+cost**2)*(1-2*sinp**2)
                  ftt2 = (2*sinp**2-1)
                  ftp = cost*2*sinp*cosp
                  fpp1 = -ftt1
                  fpp2 = cost**2*(1-2*sinp**2)
                  drt = 4*sint*(1-2*sinp**2)
                  drp = -4*sint*cost*2*sinp*cosp
                  dtt = -2*cost*(1-2*sinp**2)
                  dtp = (2-sint**2)*2*sinp*cosp
                  dpp = 2*cost*(1-2*sinp**2)
               else
                  write(*,*) 'teukwaves: m should be one of -2,-1,0,1,2'
                  STOP
               endif
               
c     parity
               if (iparity.eq.0) then    
                  grr = 1 + ca*frr
                  grt = cb*frt*rp
                  grp = cb*frp*(xp**2+yp**2)**0.5
                  gtt = rp**2*(1 + cc*ftt1 + ca*ftt2)
                  gtp = (ca - 2*cc)*ftp*rp*(xp**2+yp**2)**0.5
                  gpp = (xp**2+yp**2)*(1 + cc*fpp1 + ca*fpp2)
               else
                  grr = 1.
                  grt = ck*drt*rp
                  grp = ck*drp*rp*sint
                  gtt = (1+cl*dtt)*rp**2
                  gtp = cl*dtp*rp**2*sint
                  gpp = (1+cl*dpp)*rp**2*sint**2
               endif

c     Finally, convert the spherical components to cartesian

c     define derivatives drx = (dr/dx)               
               drx = xp/rp
               dry = yp/rp
               drz = zp/rp
      
               dtx = xp*zp/(rp**2*sqrt(xp**2+yp**2))
               dty = yp*zp/(rp**2*sqrt(xp**2+yp**2))
               dtz = (zp**2/rp**2-1)/sqrt(xp**2+yp**2)
               
               dpx = -yp/(xp**2+yp**2)
               dpy = xp/(xp**2+yp**2)
               dpz = 0.

c     cartesian g
               gxx(i,j,k) = drx**2*grr+2.*drx*dtx*grt
     &              +2.*drx*dpx*grp+dtx**2*gtt
     &              +2.*dtx*dpx*gtp+dpx**2*gpp
               
               gyy(i,j,k) = dry**2*grr+2.*dry*dty*grt
     &              +2.*dry*dpy*grp+dty**2*gtt
     &              +2.*dty*dpy*gtp+dpy**2*gpp
 
               gzz(i,j,k) = drz**2*grr+2.*drz*dtz*grt
     &              +2.*drz*dpz*grp+dtz**2*gtt
     &              +2.*dtz*dpz*gtp+dpz**2*gpp
               
               gxy(i,j,k) = drx*dry*grr+(drx*dty+dtx*dry)*grt
     &              +(drx*dpy+dpx*dry)*grp
     &              +dtx*dty*gtt+(dtx*dpy+dpx*dty)*gtp+dpx*dpy*gpp
               
               gxz(i,j,k) = drx*drz*grr+(drx*dtz+dtx*drz)*grt
     &              +(drx*dpz+dpx*drz)*grp
     &              +dtx*dtz*gtt+(dtx*dpz+dpx*dtz)*gtp+dpx*dpz*gpp
               
               gyz(i,j,k) = dry*drz*grr+(dry*dtz+dty*drz)*grt
     &              +(dry*dpz+dpy*drz)*grp
     &              +dty*dtz*gtt+(dty*dpz+dpy*dtz)*gtp+dpy*dpz*gpp
                


c     for the ext. curv. well need...
               cadot = 3*( (fpc-fnc)/rp**3 -3*(fnb+fpb)/rp**4 
     &              +3*(fpa-fna)/rp**5 )
               cbdot = -( -(fnd+fpd)/rp**2 +3*(fpc-fnc)/rp**3 
     &              -6*(fnb+fpb)/rp**4
     &              +6*(fpa-fna)/rp**5 ) 
               ccdot = ( (fpe-fne)/rp -2*(fnd+fpd)/rp**2 +9*(fpc-fnc)/rp**3
     &              -21*(fnb+fpb)/rp**4 +21*(fpa-fna)/rp**5 )/4
               ckdot = (fpc-fnc)/rp**2-3*(fpb+fnb)/rp**3+3*(fpa-fna)/rp**4
               cldot = -(fpd+fnd)/rp+2*(fpc-fnc)/rp**2-3*(fpb+fnb)/rp**3+
     &              3*(fpa-fna)/rp**4

c     parity
               if (iparity.eq.0) then    
                  hrr = -0.5*(cadot*frr)
                  hrt = -0.5*(cbdot*frt*rp)
                  hrp = -0.5*(cbdot*frp*(xp**2+yp**2)**0.5)
                  htt = -0.5*(rp**2*(ccdot*ftt1 + cadot*ftt2))
                  htp = -0.5*((cadot - 2*ccdot)*ftp*rp*(xp**2+yp**2)**0.5)
                  hpp = -0.5*((xp**2+yp**2)*(ccdot*fpp1 + cadot*fpp2))
               else
                  hrr = 0.
                  hrt = -0.5*(ckdot*drt*rp)
                  hrp = -0.5*(ckdot*drp*rp*sint)
                  htt = -0.5*((cldot*dtt)*rp**2)
                  htp = -0.5*(cldot*dtp*rp**2*sint)
                  hpp = -0.5*((cldot*dpp)*rp**2*sint**2)
               endif
               
c     time symmetry
               hxx(i,j,k) = drx**2*hrr+2.*drx*dtx*hrt
     &              +2.*drx*dpx*hrp+dtx**2*htt
     &              +2.*dtx*dpx*htp+dpx**2*hpp

               hyy(i,j,k) = dry**2*hrr+2.*dry*dty*hrt
     &              +2.*dry*dpy*hrp+dty**2*htt
     &              +2.*dty*dpy*htp+dpy**2*hpp
               
               hzz(i,j,k) = drz**2*hrr+2.*drz*dtz*hrt
     &              +2.*drz*dpz*hrp+dtz**2*htt
     &              +2.*dtz*dpz*htp+dpz**2*hpp
               
               hxy(i,j,k) = drx*dry*hrr+(drx*dty+dtx*dry)*hrt
     &              +(drx*dpy+dpx*dry)*hrp
     &              +dtx*dty*htt+(dtx*dpy+dpx*dty)*htp+dpx*dpy*hpp
               
               hxz(i,j,k) = drx*drz*hrr+(drx*dtz+dtx*drz)*hrt
     &              +(drx*dpz+dpx*drz)*hrp
     &              +dtx*dtz*htt+(dtx*dpz+dpx*dtz)*htp+dpx*dpz*hpp
               
               hyz(i,j,k) = dry*drz*hrr+(dry*dtz+dty*drz)*hrt
     &              +(dry*dpz+dpy*drz)*hrp
     &              +dty*dtz*htt+(dty*dpz+dpy*dtz)*htp+dpy*dpz*hpp
                

               
c$$$               hxx(i,j,k) = 0.
c$$$               hxy(i,j,k) = 0.
c$$$               hxz(i,j,k) = 0.
c$$$               hyy(i,j,k) = 0.
c$$$               hyz(i,j,k) = 0.
c$$$               hzz(i,j,k) = 0.


            enddo
         enddo
      enddo

 
      return
      end
      
