c/*@@
c  @file      IDAxiBrillBH.F
c  @date      
c  @author    
c  @desc 
c  
c  @enddesc 
c@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

c/*@@
c  @routine    IDAxiBrillBH
c  @date       
c  @author     
c  @desc 
c  
c  @enddesc 
c  @calls     
c  @calledby   
c  @history 
c
c  @endhistory 
c@@*/
      
c     Need include file from Einstein
#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine IDAxiBrillBH(CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

c     Perhaps this and others should go into cctk.h
      integer CCTK_Equals
      
      
      real*8 axibheps, rmax, dq, deta
      integer levels,id5,id9,idi,idg,ier
      real*8, allocatable :: cc(:,:),ce(:,:),cw(:,:),cn(:,:),cs(:,:),
     $     rhs(:,:),psi2d(:,:),detapsi2d(:,:),dqpsi2d(:,:),
     $     detaetapsi2d(:,:),detaqpsi2d(:,:),dqqpsi2d(:,:)
      real*8, allocatable :: etagrd(:),qgrd(:)
      real*8, allocatable :: eta(:,:,:),abseta(:,:,:),sign_eta(:,:,:),
     $     q(:,:,:),phi(:,:,:)
      real*8, allocatable :: psi2dv(:,:,:),detapsi2dv(:,:,:),
     $     dqpsi2dv(:,:,:),detaetapsi2dv(:,:,:),detaqpsi2dv(:,:,:),
     $     dqqpsi2dv(:,:,:)
      parameter(axibheps = 1.0d-12)
      real*8 ep1,ep2
      real*8 o1,o2,o3,o4,o5,o6,o7,o8,o9
      real*8 o10,o11,o12,o13,o14,o15,o16,o17,o18,o19
      real*8 o20,o21,o22,o23,o24,o25,o26,o27,o28,o29
      real*8 o30,o31,o32,o33,o34,o35,o36,o37,o38,o39
      real*8 o40,o41,o42,o43,o44,o45,o46,o47,o48,o49
      real*8 o50,o51,o52,o53,o54,o55,o56,o57,o58,o59
      real*8 o60,o61,o62,o63,o64,o65,o66,o67,o68,o69
      real*8 o70,o71,o72,o73,o74,o75,o76,o77,o78,o79
      real*8 o80,o81,o82,o83,o84,o85,o86,o87,o88,o89
      real*8 o90,o91,o92,o93,o94,o95,o96,o97,o98,o99
      integer i22
      real*8 pi
      real*8 adm
      integer :: nx,ny,nz
      integer i,j,k,nquads
      integer npoints,handle,ierror
      integer neb, nqb
      
      pi = 4.0d0*atan(1.0d0)
      
c     Set up the grid spacings
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)
      
c     Brill wave parameters

      print *,"Brill wave + BH Axisymmetric solve"
      write (*,123)amp,eta0,sigma,n
      print *,'etamax=',etamax
 123  format(1x,'Pars : Amp',f8.5,'  eta0',f8.5,'  sigma',f8.5,'  n ',i3)

c     Solve on this sized cartesian grid
c     2D grid size NE x NQ
c     Add 2 zones for boundaries...
c     21/11/00 TR: dont change parameters in place
c                  but keep a copy in local variables
c                  Otherwise the changed parameters cause trouble
c                  after recovery.
      neb = ne+2
      nqb = nq+2
      ! do I need to call free?
      allocate(cc(neb,nqb),ce(neb,nqb),cw(neb,nqb),cn(neb,nqb),cs(neb,nqb),
     $     rhs(neb,nqb),psi2d(neb,nqb),detapsi2d(neb,nqb),dqpsi2d(neb,nqb),
     $     detaetapsi2d(neb,nqb),detaqpsi2d(neb,nqb),dqqpsi2d(neb,nqb),
     $     etagrd(neb),qgrd(nqb))
      allocate(eta(nx,ny,nz),abseta(nx,ny,nz),sign_eta(nx,ny,nz),
     $     q(nx,ny,nz),phi(nx,ny,nz),
     $     psi2dv(nx,ny,nz),detapsi2dv(nx,ny,nz),dqpsi2dv(nx,ny,nz),
     $     detaetapsi2dv(nx,ny,nz),detaqpsi2dv(nx,ny,nz),
     $     dqqpsi2dv(nx,ny,nz))
c Initialize some arrays
      psi2d = 1.0d0
      detapsi2d = 0.0d0

      nquads = 2
      dq = nquads*0.5*pi/(nqb-2)
      deta = etamax/(neb-3)

      do j=1,nqb
	qgrd(j) = (j-1.5)*dq
	do i=1,neb
	  etagrd(i) = (i-2)*deta
#include "CactusEinstein/IDAxiBrillBH/src/bhbrill.x"
	enddo
      enddo
c     Boundary conditions
      do j=1,nqb
	ce(2,j)=ce(2,j)+cw(2,j)
	cw(2,j)=0.0

	cw(neb-1,j)=cw(neb-1,j)+ce(neb-1,j)
	cc(neb-1,j)=cc(neb-1,j)-deta*ce(neb-1,j)
	ce(neb-1,j)=0.0

      enddo
      do i=1,neb
	cc(i,2)=cc(i,2)+cs(i,2)
	cs(i,2)=0.0
	cc(i,nqb-1)=cc(i,nqb-1)+cn(i,nqb-1)
	cn(i,nqb-1)=0.0
      enddo

c     Do the solve
      print *, "     Calling axisymmetric solver"
      call mgparm (levels,5,id5,id9,idi,idg,neb,nqb)
      call mg5 (neb,2,neb-1,nqb,2,nqb-1,
     $          cc,cn,cs,cw,ce,psi2d,rhs,
     $          id5,id9,idi,idg,1,axibheps,rmax,ier)
      print *, "     Solve complete"
c     The solution is now available.
c     Debugging is needed, a stop statement should
c     be called if the IVP solve is not successful

      if(ier .ne. 0) stop 'bad solution to brill wave problem'
      print *,'rmax = ',rmax
      print *,'axibheps = ',axibheps
      print *,'psi2d = ',maxval(psi2d),' ',minval(psi2d)

      ep2 = 0.0
      do j=2,nqb-1
        do i=2,neb-1
          ep1 = rhs(i,j)-psi2d(i,j)*cc(i,j)-psi2d(i,j+1)*cn(i,j)-psi2d(i,j-1)*cs(i,j)-
     &      psi2d(i+1,j)*ce(i,j)-psi2d(i-1,j)*cw(i,j)
          ep2 = max(abs(ep1),ep2)
        enddo
      enddo
      print *,'Resulting eps =',ep2

      ! what a pain all this is....
      do j=1,nqb
         psi2d(1,j)=psi2d(3,j)
         psi2d(neb,j)=-deta*psi2d(neb-1,j)+psi2d(neb-2,j)
      enddo
      do i=1,neb
        psi2d(i,1)=psi2d(i,2)
        psi2d(i,nqb)=psi2d(i,nqb-1)
      enddo
c      goto 111
      do j=2,nqb-1
         do i=2,neb-1
            dqpsi2d(i,j)=0.5*(psi2d(i,j+1)-psi2d(i,j-1))/dq
            dqqpsi2d(i,j)=(psi2d(i,j+1)+psi2d(i,j-1)-2.*psi2d(i,j))/dq**2
            detapsi2d(i,j)=sinh(0.5*etagrd(i))+0.5*(psi2d(i+1,j)-psi2d(i-1,j))/deta
            detaetapsi2d(i,j)=0.5*cosh(0.5*etagrd(i))+
     $           (psi2d(i+1,j)+psi2d(i-1,j)-2.*psi2d(i,j))/deta**2
         enddo
      enddo
      do j=1,nqb
        detapsi2d(1,j)=-detapsi2d(3,j)
        detapsi2d(neb,j)=detapsi2d(neb-2,j) ! simplified

        detaetapsi2d(1,j)=detaetapsi2d(3,j)
        detaetapsi2d(neb,j)=detaetapsi2d(neb-2,j) ! simplified...

        dqqpsi2d(1,j)=dqqpsi2d(3,j)
        dqqpsi2d(neb,j)=dqqpsi2d(neb-2,j) ! simplified

        dqpsi2d(1,j)=dqpsi2d(3,j)
        dqpsi2d(neb,j)=-dq*dqpsi2d(neb-1,j)+dqpsi2d(neb-2,j)
      enddo
      do i=1,neb
        detapsi2d(i,1)=detapsi2d(i,2)
        detapsi2d(i,nqb)=detapsi2d(i,nqb-1)

        detaetapsi2d(i,1)=detaetapsi2d(i,2)
        detaetapsi2d(i,nqb)=detaetapsi2d(i,nqb-1)

        dqqpsi2d(i,1)=dqqpsi2d(i,2)
        dqqpsi2d(i,nqb)=dqqpsi2d(i,nqb-1)

        dqpsi2d(i,1)=-dqpsi2d(i,2)
        dqpsi2d(i,nqb)=-dqpsi2d(i,nqb-1)
      enddo
      do j=2,nqb-1
        do i=2,neb-1
          detaqpsi2d(i,j)=0.5*(detapsi2d(i,j+1)-detapsi2d(i,j-1))/dq
        enddo
      enddo
      do j=1,nqb
        detaqpsi2d(1,j)=-detaqpsi2d(3,j)
        detaqpsi2d(neb,j)=detaqpsi2d(neb-2,j) ! simplified
      enddo
      do i=1,ne
        detaqpsi2d(i,1)=-detaqpsi2d(i,2)
        detaqpsi2d(i,nqb)=-detaqpsi2d(i,nqb-1)
      enddo
      do j=1,nqb
        psi2d(:,j)=psi2d(:,j)+2.*cosh(0.5*etagrd)
      enddo

c     Now evaluate each of the following at x(i,j,k), y(i,j,k) and 
c     z(i,j,k) where i,j,k go from 1 to nx,ny,nz
      
c     Conformal factor

      eta = 0.5d0 * dlog (x**2 + y**2 + z**2)
      abseta = abs (eta)
      q = datan2 (sqrt (x**2 + y**2), z)
      phi = datan2 (y, x)

      do k=1,nz
         do j=1,ny
            do i=1,nx
c               eta(i,j,k) = 0.5d0*dlog(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)
c               abseta(i,j,k) = abs(eta(i,j,k))
               if(eta(i,j,k) .lt. 0)then
                  sign_eta(i,j,k) = -1
               else
                  sign_eta(i,j,k) = 1
               endif
c               q(i,j,k) = datan2(sqrt(x(i,j,k)**2+y(i,j,k)**2),z(i,j,k))
c                                    TYPO HERE ???????????
c                                                 |
c                                                 |
c               phi(i,j,k)= datan2(y(i,j,k),x(i,j,k))
            enddo
         enddo
      enddo

      npoints = nx*ny*nz

!     Interpolator handle.

      handle = -1

      if (interpolation_order .eq. 1) then
        call CCTK_InterpHandle (handle, "first-order uniform cartesian")
      else if (interpolation_order .eq. 2) then
        call CCTK_InterpHandle (handle, "second-order uniform cartesian")
      else if (interpolation_order .eq. 3) then
        call CCTK_InterpHandle (handle, "third-order uniform cartesian")
      endif

      if (handle .lt. 0) then
        call CCTK_WARN (0, "Couldn't get handle for interpolation operator")
      endif

      call CCTK_InterpLocal (ierror, cctkGH, handle, npoints, 2, 6, 6,
     $                       neb, nqb, etagrd, qgrd,
     $                       CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
     $                       abseta, q,
     $                       CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $                       psi2d, detapsi2d, dqpsi2d, detaetapsi2d,
     $                       detaqpsi2d, dqqpsi2d,
     $                       CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $                       CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $                       CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $                       psi2dv, detapsi2dv, dqpsi2dv, detaetapsi2dv,
     $                       detaqpsi2dv, dqqpsi2dv,
     $                       CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $                       CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $                       CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL)

      psi = psi2dv * exp (-0.5 * eta)
      detapsi2dv = sign_eta * detapsi2dv
      detaqpsi2dv = sign_eta * detaqpsi2dv

      do k=1,nz
         do j=1,ny
            do i=1,nx
               
c               psi(i,j,k) = psi2dv(i,j,k)*exp(-0.5*eta(i,j,k))
c               detapsi2dv(i,j,k) = sign_eta(i,j,k)*detapsi2dv(i,j,k)
               
c     psix = \partial psi / \partial x / psi
#include "CactusEinstein/IDAxiBrillBH/src/psi_1st_deriv.x"
               
c               detaqpsi2dv(i,j,k) = sign_eta(i,j,k)*detaqpsi2dv(i,j,k)
               
c     psixx = \partial^2\psi / \partial x^2 / psi
#include "CactusEinstein/IDAxiBrillBH/src/psi_2nd_deriv.x"
            enddo
         enddo
      enddo
      
c     Conformal metric
c     gxx = ...
      
c     Derivatives of the metric
c     dxgxx = 1/2 \partial gxx / \partial x

      do k=1,nz
         do j=1,ny
            do i=1,nx
c THESE WERE ALREADY CALCULATED ABOVE !!!
c               eta(i,j,k) = 0.5d0*dlog(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2) 
c               q(i,j,k) = datan2(sqrt(x(i,j,k)**2+y(i,j,k)**2),z(i,j,k))
c               phi(i,j,k) = datan2(y(i,j,k),x(i,j,k))
#include "CactusEinstein/IDAxiBrillBH/src/gij.x"
            enddo
         enddo
      enddo

c     Curvature
      kxx = 0.0D0
      kxy = 0.0D0
      kxz = 0.0D0
      kyy = 0.0D0
      kyz = 0.0D0
      kzz = 0.0D0
      
 111  continue
c     Set ADM mass
      i = neb-15
      adm = 0.0
      do j=2,nqb-1
        adm=adm+(psi2d(i,j)-(psi2d(i+1,j)-psi2d(i-1,j))/deta)*exp(0.5*
     $        etagrd(i))
      enddo
      adm=adm/(nqb-2)
      print *,'ADM mass: ',adm
      if (CCTK_Equals(initial_lapse,"schwarz")==1) then
         write (*,*)"Initial with schwarzschild-like lapse"
         write (*,*)"using alp = (2.*r - adm)/(2.*r+adm)."
         alp = (2.*r - adm)/(2.*r+adm)
      endif
      
      conformal_state = CONFORMAL_METRIC
      
      deallocate(cc,ce,cw,cn,cs,rhs,psi2d,detapsi2d,dqpsi2d,
     $     detaetapsi2d,detaqpsi2d,dqqpsi2d,
     $     etagrd,qgrd,
     $     eta,abseta,sign_eta,q,phi,psi2dv,detapsi2dv,dqpsi2dv,
     $     detaetapsi2dv,detaqpsi2dv,dqqpsi2dv)

      return
      end

