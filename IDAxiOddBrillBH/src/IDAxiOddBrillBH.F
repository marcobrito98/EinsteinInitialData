c/*@@
c  @file      IDAxiOddBrillBH.F
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
c  @routine    IDAxiOddBrillBH
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

      subroutine IDAxiOddBrillBH(CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

c     Perhaps this and others should go into cctk.h
      integer CCTK_Equals   

      real*8 deta,dq
      real*8, allocatable :: ac(:,:),ae(:,:),aw(:,:),an(:,:),as(:,:),
     $     rhs(:,:),ksq(:,:), psi2dv(:,:),dpsi2dv(:,:),ddpsi2dv(:,:),
     $     detapsisph(:,:),dqpsisph(:,:),detaetapsisph(:,:),
     $     detaqpsisph(:,:),dqqpsisph(:,:),
     $     etagrd(:),qgrd(:)
      real*8 o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,
     $     o11,o12,o13,o14,o15,o16,o17,o18,o19,
     $     o20,o21,o22,o23,o24,o25,o26,o27,o28,o29,
     $     o30,o31,o32,o33,o34,o35,o36,o37,o38,o39,
     $     o40,o41,o42,o43,o44,o45,o46,o47,o48,o49,
     $     o50,o51
      real*8 t1,t2,t3,t4
      real*8 gtil,dngtil,dnngtil,dnnngtil,dnnnngtil,dnnnnngtil,
     $     dnnnnnngtil,dnnnnnnngtil
      real*8 rhsmax,rmax,odd_get2d,adm,Jmom,mass,a,R,rBL
      real*8,parameter :: rbh_tol = 1.0d-7,rbh_eps = 1.0d-10
      integer,parameter :: itmax = 100 
      real*8 pi,zero,one,two
      CCTK_REAL :: xmin,xmax,ymin,ymax,zmin,zmax
      integer :: ne, nq
      integer :: nx,ny,nz
      integer :: i,j,k,it,ier,nquads,ntries
      integer :: npoints,handle,ierror

      pi = 4.0d0*atan(1.0d0)
c     DO NOT use integer*4
      ne = neta
      nq = ntheta

c     Set up the grid spacings
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)      

c     Get parameters
c     --------------
      if (brandt_seidel == 1) then
         write(*,*) '"Solve Bradt-Seidel Initial Data Sets"'
      endif

      write(*,*) 'Brill Wave: amplitude :',amp
      write(*,*) 'Brill Wave: center :',eta0
      write(*,*) 'Brill Wave: sigma :',sigma
      write(*,*) 'Bowen & York Momenta:  :',byJ
      write(*,*) 'Brill Wave: sin^n theta :',n

      conformal_state = CONFORMAL_METRIC
c     Sovle on this sized cartesian grid
c     3D grid size NE x NT x NP
c     Add 2 zones for eta coordinate and 4 for theta
      ne = ne+2
      nq = nq+2
c     
      allocate(ac(ne,nq),ae(ne,nq),aw(ne,nq),an(ne,nq),as(ne,nq),
     $     rhs(ne,nq),ksq(ne,nq),psi2dv(ne,nq),dpsi2dv(ne,nq),
     $     ddpsi2dv(ne,nq),detapsisph(ne,nq),dqpsisph(ne,nq),
     $     detaetapsisph(ne,nq),detaqpsisph(ne,nq),dqqpsisph(ne,nq),
     $     etagrd(ne),qgrd(nq))
c     
c     Initialize some array
c     
      nquads = 2      
      dq = nquads*0.5*pi/(nq-2)
      deta = etamax/(ne-3)
      
      do j = 1,nq
         qgrd(j) = (j-1.5)*dq
      enddo
      do i=1,ne
         etagrd(i) = (i-2)*deta
      enddo
      
c     
c     Specify the initial guess for the conformal factor:
c     
      dpsi2dv = 0.
      ddpsi2dv = 0.
      do j = 1,nq
         do i = 1,ne
            psi2dv(i,j) = 2.*cosh(0.5*etagrd(i))
         enddo
      enddo
c     
c     Compute K^2:
c     
      do j = 1,nq
         do i = 1,ne
#include "Development/IDAxiOddBrillBH/src/gauss.x"
            if(brandt_seidel == 1) then
#include "Development/IDAxiOddBrillBH/src/ksq_bs.x"
	    else if (sergio == 1) then
#include "Development/IDAxiOddBrillBH/src/kerr.x"
#include "Development/IDAxiOddBrillBH/src/ksq_sergio.x"
            else
#include "Development/IDAxiOddBrillBH/src/ksq_axi_rdbh.x"
            endif
         enddo
      enddo
c     
c     Solve the Hamiltonian constraint for the conformal factor:
c     

      ntries = 0

      do it = 1,itmax+1
         
         ac = 0.
         ae = 0.
         aw = 0.
         an = 0.
         as = 0.
         rhs = 0.

c     
         do j = 2,nq-1
            do i = 2,ne-1
               rhs(i,j) =
     $              (dpsi2dv(i+1,j)-2.*dpsi2dv(i,j)+dpsi2dv(i-1,j))/
     $              deta**2+(dpsi2dv(i,j+1)-2.*dpsi2dv(i,j)+
     $              dpsi2dv(i,j-1))/dq**2+0.5*(dpsi2dv(i,j+1)-
     $              dpsi2dv(i,j-1))/(dq*tan(qgrd(j)))-0.25*dpsi2dv(i,j)+
     $              0.125*ksq(i,j)/(psi2dv(i,j)+dpsi2dv(i,j))**7
            enddo
         enddo
c     
         rhs = -rhs
c     
c     See if the residual is small enough:
c     
         rhsmax = 0.
c     
         do j = 2,nq-1
            do i = 2,ne-1
               rhsmax = max(rhsmax,abs(rhs(i,j)))
            enddo
         enddo

         if (verbose == 1) then
            ntries = ntries+1
            print*,'-----------Retry', ntries
            print*,'residual =', rhsmax
            print*,'--------------------------------------------------'
         endif

         if (rhsmax.lt.rbh_tol) goto 110
c     
c     Compute stencil coefficients:
c     
         do j = 2,nq-1
            do i = 2,ne-1
               ac(i,j) = -2./deta**2-2./dq**2-sin(qgrd(j))**2-0.25-
     $              0.875*ksq(i,j)/(psi2dv(i,j)+ 
     $              dpsi2dv(i,j))**8
c     
               ae(i,j) = 1./deta**2
c     
               aw(i,j) = 1./deta**2
c     
               an(i,j) = 1./dq**2 + 0.5/(tan(qgrd(j))*dq)
c     
               as(i,j) = 1./dq**2 - 0.5/(tan(qgrd(j))*dq)
c     
            enddo
         enddo
         
c     
c     Apply boundary conditions:
c     
c     i=2:
c     
         do j = 3,nq-2
            ae(2,j) = ae(2,j) + aw(2,j)
            aw(2,j) = 0.
c     
c     i=ne-1:
c     
            ac(ne-1,j) = ac(ne-1,j) + 4.*ae(ne-1,j)/(3.+deta)
            aw(ne-1,j) = aw(ne-1,j) - ae(ne-1,j)/(3.+deta)
            ae(ne-1,j) = 0.
         enddo
c     
c     j=2:
c     
         do i = 3,ne-2
            ac(i,2) = ac(i,2) + as(i,2)
            as(i,2) = 0.
c     
c     j=nq-1:
c     
            ac(i,nq-1) = ac(i,nq-1) + an(i,nq-1)
            an(i,nq-1) = 0.
         enddo
         
c     
c     i=2, j=2:
c     
         ae(2,2) = ae(2,2) + aw(2,2)
         ac(2,2) = ac(2,2) + as(2,2)
         aw(2,2) = 0.
         as(2,2) = 0.
c     
c     i=2, j=nq-1:
c     
         ae(2,nq-1) = ae(2,nq-1) + aw(2,nq-1)
         ac(2,nq-1) = ac(2,nq-1) + an(2,nq-1)
         aw(2,nq-1) = 0.
         an(2,nq-1) = 0.
c     
c     i=ne-1, j=2:
c     
         ac(ne-1,2) = ac(ne-1,2) + 4.*ae(ne-1,2)/(3.+deta) + 
     $        as(ne-1,2)
         aw(ne-1,2) = aw(ne-1,2) - ae(ne-1,2)/(3.+deta)
         ae(ne-1,2) = 0.
         as(ne-1,2) = 0.
c     
c     i=ne-1, j=nq-1:
c     
         ac(ne-1,nq-1) = ac(ne-1,nq-1) + 4.*ae(ne-1,nq-1)/(3.+deta) +
     $        an(ne-1,nq-1)
         aw(ne-1,nq-1) = aw(ne-1,nq-1) - ae(ne-1,nq-1)/(3.+deta)
         ae(ne-1,nq-1) = 0.
         an(ne-1,nq-1) = 0.
c     
         call bicgst2d(ac,ae,aw,an,as,ddpsi2dv,rhs,rbh_eps
     $        ,rmax,ier,ne,nq)
         if (rmax.gt.1.0e-10) then
            write(*,*) '***WARNING: bicgst3d did not converge.'
         endif
         if (ier.eq.-1) then
            write(*,*) '***WARNING: ier=-1'
         endif
         
c     
c     Now, apply boundary conditions to ddpsi2dv:
c     
         do j = 1,nq
            ddpsi2dv(1,j) = ddpsi2dv(3,j)
            ddpsi2dv(ne,j) = (4.*ddpsi2dv(ne-1,j)-ddpsi2dv(ne-2,j))/(3.+deta)
         enddo
         do i = 1,ne
            ddpsi2dv(i,1) = ddpsi2dv(i,2)
            ddpsi2dv(i,nq) = ddpsi2dv(i,nq-1)
         enddo
         do j = 1,nq
            do i = 1,ne
               ddpsi2dv(i,j) = ddpsi2dv(i,j)
               ddpsi2dv(i,j) = ddpsi2dv(i,j)
            enddo
         enddo
c     
c     Update dpsi2dv:
c     
         do j = 1,nq
            do i = 1,ne
               dpsi2dv(i,j) = dpsi2dv(i,j) + ddpsi2dv(i,j)
            enddo
         enddo
         
c     
      enddo
c     
      write(*,*) 'It did not converge.'
      stop
c     
 110  continue
      write(*,*) '--------------------------------------------------'
      print*,'Converge at Residual = ',rhsmax
c     
c     Here, compute the derivatives of the spherical conformal factor
c     
      do j = 1,nq
         do i = 2,ne-1
            detapsisph(i,j)=0.5*(dpsi2dv(i+1,j)-dpsi2dv(i-1,j))/
     $           deta + sinh(0.5*etagrd(i))
         enddo
         detapsisph(1,j) = -detapsisph(3,j)
      enddo
c     
      do j = 2,nq-1
         do i = 1,ne
            dqpsisph(i,j)=0.5*(dpsi2dv(i,j+1)-dpsi2dv(i,j-1))/dq
         enddo
      enddo
      do i = 1,ne
         dqpsisph(i,1) = -dqpsisph(i,2)
         dqpsisph(i,nq) = -dqpsisph(i,nq-1)
      enddo
      
c     
      do j = 1,nq
         do i = 2,ne-1
            detaetapsisph(i,j)=(dpsi2dv(i+1,j)-2.*dpsi2dv(i,j)+dpsi2dv(i-1,j))/
     $           deta**2 + sqrt(0.25)*cosh(0.5*etagrd(i))
         enddo
         detaetapsisph(1,j) = detaetapsisph(3,j)
      enddo
      
c     
      do j = 2,nq-1
         do i = 1,ne
            detaqpsisph(i,j)=0.5*(detapsisph(i,j+1)-
     $           detapsisph(i,j-1))/dq
         enddo
      enddo
      do i = 1,ne
         detaqpsisph(i,1) = -detaqpsisph(i,2)
         detaqpsisph(i,nq) = -detaqpsisph(i,nq-1)
      enddo
c     
      do j = 2,nq-1
         do i = 1,ne
            dqqpsisph(i,j)=0.5*(dqpsisph(i,j+1)-dqpsisph(i,j-1))/dq
         enddo
      enddo
      do i = 1,ne
         dqqpsisph(i,1) = dqqpsisph(i,2)
         dqqpsisph(i,nq) = dqqpsisph(i,nq-1)
      enddo
c     
      do j = 1,nq
         psi2dv(:,j)=dpsi2dv(:,j)+2.0*cosh(0.5*etagrd)
      enddo
c     
c     Now compute on the Cartesian coordinate.
c     
c     Compute eta,q,phi at the each points of cartesian grid
         
      eta = 0.5d0 * dlog(x**2 + y**2 + z**2)
      abseta = abs (eta)
      q = datan2 (sqrt (x**2 + y**2),z)
      phi = datan2 (y, x)
      
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
      
               if(eta(i,j,k) .lt. 0)then
                  sign_eta(i,j,k) = -1
               else
                  sign_eta(i,j,k) = 1
               endif
            enddo
         enddo
      enddo

     
c     Find the local origin of the spatial coordinates
c     ------------------------------------------------

      npoints = nx*ny*nz

c     Interpolator handle.

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
      
      call CCTK_InterpLocal  (ierror,cctkGH,handle,npoints,2,6,6,
     $     ne,nq,etagrd,qgrd,
     $     CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $     abseta,q,
     $     CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $     psi2dv,detapsisph,dqpsisph,detaetapsisph,detaqpsisph,
     $     dqqpsisph,
     $     CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $     CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $     psi2d,detapsi2d,dqpsi2d,detaetapsi2d,detaqpsi2d,
     $     dqqpsi2d,
     $     CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,
     $     CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL,CCTK_VARIABLE_REAL)
      
      psi = psi2d*exp(-0.5*eta)
      detapsi2d = sign_eta * detapsi2d
      detaqpsi2d = sign_eta * detaqpsi2d
  
      do k=1,nz
         do j=1,ny
            do i=1,nx
               
c     psix = \partial psi / \partial x / psi
#include "Development/IDAxiOddBrillBH/src/psi_1st_deriv.x"
               
c     psixx = \partial^2\psi / \partial x^2 / psi
#include "Development/IDAxiOddBrillBH/src/psi_2nd_deriv.x"

            enddo
         enddo
      enddo
      
c     Conformal metric (here, we define conformally flat)
c     gxx = 1...
      
c     Derivatives of the metric
c     dxgxx = 1/2 \partial gxx / \partial x = 0...
c     
      do k= 1,nz
         do j= 1,ny
            do i= 1,nx
#include "Development/IDAxiOddBrillBH/src/gij.x"
            enddo
         enddo
      enddo
c     Extrinsic Curvture:Cactus only needs physical extrinsic curvature
      do k= 1,nz
         do j= 1,ny
            do i= 1,nx 
#include "Development/IDAxiOddBrillBH/src/cgauss.x"
               if (brandt_seidel == 1) then
#include "Development/IDAxiOddBrillBH/src/kij_bs.x"
	       else if (sergio == 1) then
#include "Development/IDAxiOddBrillBH/src/ckerr.x"
#include "Development/IDAxiOddBrillBH/src/kij_sergio.x"
               else
#include "Development/IDAxiOddBrillBH/src/kij_axi_rdbh.x"
               endif
            enddo
         enddo
      enddo
      
      kxx = kxx/psi**2
      kxy = kxy/psi**2
      kxz = kxz/psi**2
      kyy = kyy/psi**2
      kyz = kyz/psi**2
      kzz = kzz/psi**2

c     Set ADM mass
      i = ne-15
      adm = 0.0
      do j=2,nq-1
         adm=adm+(psi2dv(i,j)-(psi2dv(i+1,j)-psi2dv(i-1,j))/
     $        deta)*exp(0.5*etagrd(i))
      enddo
      
      adm=adm/(nq-2)
      print *,'ADM mass: ',adm

c     Set Angular momentum
      print*,'Angular momentum parameter: a/m = ', byJ/(adm**2)
c      i = ne-15
c      Jmom = 0.0
c      do j=2,nq-1
c         Jmom=Jmom+psi2dv(i,j)**6*exc31(i,j)*sin(qgrd(j))**3*dq*.5d0
c      enddo
      
c      Jmom=Jmom/(nq-2)**9
c      print *,'Angular momentum: ',Jmom
      
      if (CCTK_Equals(initial_lapse,"schwarz")==1) then 
         write (*,*)"Initial with schwarzschild-like lapse"
         write (*,*)"using alp = (2.*r - adm)/(2.*r+adm)."
         alp = (2.*r - adm)/(2.*r+adm)
      endif
      
      deallocate(ac,ae,aw,an,as,rhs,ksq,psi2dv,dpsi2dv,ddpsi2dv,
     $     detapsisph,dqpsisph,detaetapsisph,detaqpsisph,dqqpsisph,
     $     etagrd,qgrd)
      
      return
      end
      
