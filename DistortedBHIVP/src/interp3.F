c /*@@
c   @routine    interp3
c   @date       Fri Feb 14 08:46:53 1997
c   @author     Ryoji Takahashi
c   @desc 
c      Interpolates from 3D data var with coordinates x, y, and z and
c      sizes nx , ny, and nz onto 1D data out with position outx, outy
c      ,and outz nout points.
c      <p>
c      This has linear, quadratic and cubic interpolators in it.
c      Or will one day very soon.
c   @enddesc 
c   @calls     
c   @calledby   numerical_nonaxi
c@@*/

      subroutine interp3(dat,x,y,z,nx,ny,nz,out,outx,outy,outz,nout
     $     ,order)
      implicit none 
      integer nx,ny,nz,nout
      real*8 dat(nx,ny,nz), x(nx), y(ny), z(nz)
      real*8 out(nout),outx(nout),outy(nout),outz(nout)
      integer order
c Interpolation goes from ibelow to ibelow+[1,2,3] depending on order
      integer i,j,k,l,m,ibelow,jbelow,kbelow,pt
      real*8 xsym,ysym,zsym,findx,findy,findz,frac
      real*8 ydir(order+1)
      real*8 zdir(order+1)
      real*8 f(4), fp(4,4), fl(4)
      real*8 dbh_linear, dbh_quad, dbh_cubic
      real*8 dx, dy, dz, PI
      integer twobhjsad

      pi = 4.0d0*atan(1.0d0)
      
      dx = x(2) - x(1)
      dy = y(2) - y(1)
      dz = z(2) - z(1)
      
c Loop over all out points
      do pt=1,nout
         zsym = 1.0D0
         ysym = 1.0D0
         xsym = 1.0D0
c Check bounds
         findx = outx(pt)
         if (findx .lt. x(1)) then 
            write (*,*) "Below inner bound at ",pt,outx(pt)
            STOP
         endif
         if (findx .gt. x(nx)) then
            write (*,*) "Above x bounds at ",pt,outx(pt),x(nx)
            STOP
         endif
         findy = outy(pt)
         if (findy .lt. y(1)) then 
            write (*,*) "Below  y inner bound at ",pt,outy(pt)
            STOP
         endif
         if (findy .gt. y(ny)) then 
            write (*,*) "Below  y inner bound at ",pt,outy(pt)
            STOP
         endif
         findz = outz(pt)+pi
         if (findz .lt. z(1)) then 
            write (*,*) "Below  z inner bound at ",pt,outz(pt)
            STOP
         endif
         if (findz .gt. z(nz)) then 
            write (*,*) "Below  z inner bound at ",pt,outz(pt)
            STOP
         endif

c Locate ourselves in i,j space
c         do i=1,nx
c            if (x(i) .lt. findx) then
c               ibelow = i
c            endif
c         enddo
c     Assume a regular grid
         ibelow = (findx-x(1))/dx+1

         if (order .eq. 3 .and. ibelow .gt. 1) then
            ibelow = ibelow - 1
         endif
         if (ibelow + order .gt. nx) then
            ibelow = nx - order
         endif

c         do i=1,ny
c            if (y(i) .lt. findy) then
c     jbelow = i
c     endif
c         enddo
c     Assume a regular grid
         jbelow = (findy-y(1))/dy+1
         if (order .eq. 3 .and. jbelow .gt. 1) then
            jbelow = jbelow - 1
         endif
         if (jbelow + order .gt. ny) then
            jbelow = ny - order
         endif

c         do i=1,nz
c            if (z(i) .lt. findz) then
c               kbelow = i
c            endif
c         enddo
c     Assume a regular grid
         kbelow = (findz-z(1))/dz+1
         if (order .eq. 3 .and. kbelow .gt. 1) then
            kbelow = kbelow - 1
         endif
         if (kbelow + order .gt. nz) then
            kbelow = nz - order
         endif
         
c     write (*,*) "PT :",findx,findy
c     write (*,*) "SYM:",sym
c     write (*,*) "BOUND X ",ibelow,x(ibelow),x(ibelow+1)
c     write (*,*) "BOUND Y ",jbelow,y(jbelow),y(jbelow+1)
c     write (*,*) "BOUND z ",kbelow,z(kbelow),z(kbelow+1)
         
c     So do the interpolation
         if (order .eq. 1) then 
c     Interp in the x direction
!            frac = (findx-x(ibelow))/(x(ibelow+1)-x(ibelow))
            do l = 1,2
               do m = 1,2
                  f(1) = dat(ibelow,jbelow+l-1,kbelow+m-1)
                  f(2) = dat(ibelow+1,jbelow+l-1,kbelow+m-1)
                  fp(l,m) = dbh_linear(f,x(ibelow),dx,findx)
               enddo
            enddo
c     Now take our 2x2 plane and interp to the center of both
c     in the y direction
!            frac = (findy-y(jbelow))/(y(jbelow+1)-y(jbelow))
            do m = 1,2
               f(1) = fp(1,m)
               f(2) = fp(2,m)
               fl(m) = dbh_linear(f,y(jbelow),dy,findy)
            enddo
c     And finally, interp in the z direction
            out(pt) = xsym * ysym * zsym * 
     $           dbh_linear(fl,z(kbelow),dz,findz)
         else if (order .eq. 2) then
c     Load up for calls to poly2inter
            do l=1,3
               do m=1,3
                  f(1) = dat(ibelow,jbelow+l-1,kbelow+m-1)
                  f(2) = dat(ibelow+1,jbelow+l-1,kbelow+m-1)
                  f(3) = dat(ibelow+2,jbelow+l-1,kbelow+m-1)
                  fp(l,m) = dbh_quad(f,x(ibelow),dx,findx)
               enddo
            enddo
c     Now take our 2x2 plane and interp to the center of both
c     in the y direction
            do m=1,3
               f(1) = fp(1,m)
               f(2) = fp(2,m)
               f(3) = fp(3,m)
               fl(m) = dbh_quad(f,y(jbelow),dy,findy)
            enddo
c     And finally, interp in the z direction
            out(pt) = xsym * ysym * zsym *
     $           dbh_quad(fl,z(kbelow),dz,findz)
         else if (order .eq. 3) then
c     Load up for calls to cubic
            do l=1,4
               do m=1,4
                  f(1) = dat(ibelow,jbelow+l-1,kbelow+m-1)
                  f(2) = dat(ibelow+1,jbelow+l-1,kbelow+m-1)
                  f(3) = dat(ibelow+2,jbelow+l-1,kbelow+m-1)
                  f(4) = dat(ibelow+3,jbelow+l-1,kbelow+m-1)
                  fp(l,m) = dbh_cubic(f,x(ibelow),dx,findx)
               enddo
            enddo
c     Now take our 2x2 plane and interp to the center of both
c     in the y direction
            do m=1,4
               f(1) = fp(1,m)
               f(2) = fp(2,m)
               f(3) = fp(3,m)
               f(4) = fp(4,m)
               fl(m) = dbh_cubic(f,y(jbelow),dy,findz)
            enddo
c     And finally, interp in the z direction
            out(pt) = xsym * ysym * zsym * 
     $           dbh_cubic(fl,z(kbelow),dz,findz)
         else
            write (*,*) "ORDER set wrong in interp3d",order
            stop
         endif

      enddo
      return
      end

      real*8 function dbh_linear(f, x0, dx, findx)
      implicit none
      real*8 f(2),x0,dx,findx
      real*8 frac

      frac = (findx-x0)/dx
      dbh_linear = (frac)*f(2) + (1.0-frac)*f(1)

      return
      end

      real*8 function dbh_quad(f, x0, dx, findx)
      implicit none
      real*8 f(3),x0, dx, findx, dbh_quad
      real*8 f0,f1,f2
      real*8 a,b,c, dx2, x02, o2dx2
c Mathematica tells us
c     -  List(List(Rule(c,(2*dx**2*f0 + 3*dx*f0*x0 - 4*dx*f1*x0 
c     -        + dx*f2*x0 +  f0*x0**2 - 2*f1*x0**2 + f2*x0**2)
c     -        /(2*dx**2)), Rule(b,(-3*dx*f0 + 4*dx*f1 - dx*f2 - 2*f0*x0 
c     -        + 4*f1*x0 - 2*f2*x0)/(2*dx**2)),Rule(a,(f0 - 2*f1 + 
c     -        f2)/(2*dx**2))))

	f0 = f(1)
	f1 = f(2)
	f2 = f(3)
        dx2 = dx**2
        x02 = x0**2
        o2dx2 = 1.0D0/(2.0D0*dx2)

      c = (2.0D0*dx2*f0 + dx*x0*(3.0D0*f0 - 4.0D0*f1 + f2) + 
     $       x02*(f0 - 2.0D0*f1 + f2))*o2dx2
      b = (dx * (-3.0D0*f0 + 4.0D0*f1 - f2) + x0 * (- 2.0D0*f0 + 
     $     4.0D0*f1 - 2.0D0*f2))*o2dx2

      a = (f0 - 2.0D0*f1 + f2)*o2dx2

      dbh_quad = a*findx**2 + b*findx + c

      end


      real*8 function dbh_cubic(f, x0, dx, findx)
      implicit none
      real*8 a,b,c,d,dbh_cubic
      real*8 f(4),x0,dx,findx

      a = -(f(1)-3.0*f(2)+3.0*f(3)-f(4)) / (6.0*(dx**3))
      
      b = (f(1)-2.0*f(2)+f(3))/(2.0*(dx**2)) +
     $     (f(1)-3.0*f(2)+3.0*f(3)-f(4))*(dx+x0)/(2.0*(dx**3))
      
      c = ((dx**2)*(-11.0*f(1) + 18.0*f(2) -  9.0*f(3) +  2.0*f(4)) +
     $     dx*x0*  (-12.0*f(1) + 30.0*f(2) - 24.0*f(3) +  6.0*f(4)) +
     $     (x0**2)*( -3.0*f(1) +  9.0*f(2) -  9.0*f(3) +  3.0*f(4))) /
     $     (6.0*(dx**3))
      
      d = ((dx**3)*   (  6.0*f(1) ) +
     $     (dx**2)*x0*( 11.0*f(1) - 18.0*f(2) +  9.0*f(3) -  2.0*f(4)) +
     $     (x0**2)*dx*(  6.0*f(1) - 15.0*f(2) + 12.0*f(3) -  3.0*f(4)) +
     $     (x0**3)*   (  1.0*f(1) -  3.0*f(2) +  3.0*f(3) -  1.0*f(4)))/
     $     (6.0*(dx**3))
         
      dbh_cubic = ((a*findx + b)*findx + c)*findx + d
      return 
      end
