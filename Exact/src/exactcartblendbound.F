#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine exactcartblendboundary(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      logical doKij, doGij, doDs, doVs, doLapse, doShift

      integer i,j,k
      integer nx,ny,nz
      integer ninterps
      integer CCTK_Equals

      CCTK_REAL xblend, yblend, zblend
      CCTK_REAL xmin, xmax, ymin, ymax, zmin, zmax
      CCTK_REAL xfrac, yfrac, zfrac, onemxfrac, onemyfrac, onemzfrac
      CCTK_REAL oonints
      CCTK_REAL sfrac, onemsfrac
      CCTK_REAL gxxe, gyye, gzze, gxye, gyze, gxze
      CCTK_REAL kxxe, kyye, kzze, kxye, kyze, kxze
      CCTK_REAL dxgxxe, dxgyye, dxgzze, dxgxye, dxgyze, dxgxze
      CCTK_REAL dygxxe, dygyye, dygzze, dygxye, dygyze, dygxze
      CCTK_REAL dzgxxe, dzgyye, dzgzze, dzgxye, dzgyze, dzgxze
      CCTK_REAL alpe, axe, aye, aze
      CCTK_REAL betaxe,betaye,betaze
      CCTK_REAL bxxe,bxye,bxze,byxe,byye,byze,bzxe,bzye,bzze
      CCTK_REAL det, uxx, uxy, uxz, uyy, uyz, uzz
      CCTK_REAL vxe,vye,vze,sav

      CCTK_REAL :: dx,dy,dz,time,ierr

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      dx = cctk_delta_space(1)
      dy = cctk_delta_space(2)
      dz = cctk_delta_space(3)

      time = cctk_time

C     Other parameters.

      doKij = (exblend_Ks.eq.1)
      doGij = (exblend_gs.eq.1)

      doLapse  = ((exblend_gauge.eq.1).and.
     $            (CCTK_Equals("lapse_evolution_method","exact").ne.0))
      doShift  = ((exblend_gauge.eq.1).and.
     $            (CCTK_Equals(shift_evolution_method,"exact").ne.0))

      if (exblend_width.lt.0) then 
         xblend = - exblend_width*dx
         yblend = - exblend_width*dy
         zblend = - exblend_width*dz
      else
         xblend = exblend_width
         yblend = exblend_width
         zblend = exblend_width
      endif

      call CCTK_CoordRange(ierr,cctkGH,xmin,xmax,-1,"x","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,ymin,ymax,-1,"y","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,zmin,zmax,-1,"z","cart3d")

      do k=1,nz
         do j=1,ny
            do i=1,nx

c     We only do anything if in the blending region

               if (x(i,j,k) .ge. xmax - xblend .or.
     $             x(i,j,k) .le. xmin + xblend .or.
     $             y(i,j,k) .ge. ymax - yblend .or.
     $             y(i,j,k) .le. ymin + yblend .or.
     $             z(i,j,k) .ge. zmax - zblend .or.
     $             z(i,j,k) .le. zmin + zblend) then

                  call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), time,
     $                 gxxe, gyye, gzze, gxye, gyze, gxze,
     $                 kxxe, kyye, kzze, kxye, kyze, kxze,
     $                 dxgxxe, dxgyye, dxgzze, dxgxye, dxgyze, dxgxze,
     $                 dygxxe, dygyye, dygzze, dygxye, dygyze, dygxze,
     $                 dzgxxe, dzgyye, dzgzze, dzgxye, dzgyze, dzgxze,
     $                 alpe, axe, aye, aze, betaxe, betaye, betaze,
     $                 bxxe, bxye, bxze, byxe, 
     $                 byye, byze, bzxe, bzye, bzze)

c     This sucks, but we want the exact vs so we can blend them also.

                  det = -(gxze**2*gyye) 
     &                 + 2.d0*gxye*gxze*gyze 
     &                 - gxxe*gyze**2 
     &                 - gxye**2*gzze 
     &                 + gxxe*gyye*gzze
  
                  uxx=(-gyze**2 + gyye*gzze)/det
                  uxy=(gxze*gyze - gxye*gzze)/det
                  uyy=(-gxze**2 + gxxe*gzze)/det
                  uxz=(-gxze*gyye + gxye*gyze)/det
                  uyz=(gxye*gxze - gxxe*gyze)/det
                  uzz=(-gxye**2 + gxxe*gyye)/det

c     OK so 6 blending cases. If frac = 1 we get all exact

                  ninterps = 0

                  xfrac     = 0.0D0
                  onemxfrac = 0.0D0
                  yfrac     = 0.0D0
                  onemyfrac = 0.0D0
                  zfrac     = 0.0D0
                  onemzfrac = 0.0D0

                  if (x(i,j,k) .le. xmin + xblend) then 
                     xfrac = 1.0D0 - (x(i,j,k)-xmin) / xblend
                     onemxfrac = 1.0D0 - xfrac
                     ninterps = ninterps + 1
                  endif

                  if (x(i,j,k) .ge. xmax - xblend) then 
                     xfrac = 1.0D0 - (xmax - x(i,j,k)) / xblend
                     onemxfrac = 1.0D0 - xfrac
                     ninterps = ninterps + 1
                  endif

                  if (y(i,j,k) .le. ymin + yblend) then 
                     yfrac = 1.0D0 - (y(i,j,k)-ymin) / yblend
                     onemyfrac = 1.0D0 - yfrac
                     ninterps = ninterps + 1
                  endif

                  if (y(i,j,k) .ge. ymax - yblend) then 
                     yfrac = 1.0D0 - (ymax - y(i,j,k)) / yblend
                     onemyfrac = 1.0D0 - yfrac
                     ninterps = ninterps + 1
                  endif

                  if (z(i,j,k) .le. zmin + zblend) then 
                     zfrac = 1.0D0 - (z(i,j,k)-zmin) / zblend
                     onemzfrac = 1.0D0 - zfrac
                     ninterps = ninterps + 1
                  endif

                  if (z(i,j,k) .ge. zmax - zblend) then 
                     zfrac = 1.0D0 - (zmax - z(i,j,k)) / zblend
                     onemzfrac = 1.0D0 - zfrac
                     ninterps = ninterps + 1
                  endif

                  oonints = 1.0D0 / ninterps

                  if (ninterps .eq. 0 .or. ninterps .gt. 3) then 
                     print *,"NINTERPS error", ninterps
                     STOP
                  endif
                     
                  sfrac = (xfrac + yfrac + zfrac) * oonints
                  onemsfrac = 1.0D0 - sfrac

c     Once again some c-preprocessor tricks based on the whole fortran
c     space thing...

#define INTPOINT(f,v) f(i,j,k) = sfrac * v  + onemsfrac * f(i,j,k)
#define intone(f) INTPOINT(f, f e)
#define int_grp(p) \
        intone(p xx) &&\
        intone(p xy) &&\
        intone(p xz) &&\
        intone(p yy) &&\
        intone(p yz) &&\
        intone(p zz)

                  if (doGij) then 
                     int_grp(g)
                  endif

                  if (doKij) then
                     int_grp(k)
                  endif

                  if (doLapse) then 
                     intone(alp)
                  endif

                  if (doShift.and.(shift_state.ne.SHIFT_INACTIVE)) then 
                     intone(betax)
                     intone(betay)
                     intone(betaz)
                  endif

               endif            ! r > rinner

            enddo
         enddo
      enddo

      return
      end 
