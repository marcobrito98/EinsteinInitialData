c=======================================================================
c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c by Hisaaki Shinkai   shinkai@wurel.wustl.edu  19980603
c-----------------------------------------------------------------------
c  This is for maximally charged multi BH solutions such as 
c  Majumdar-Papapetrou (1947) or Kastor-Traschen (1993) solution.
c  See also doc/KTsol.tex for brief review of this solution. 
c-----------------------------------------------------------------------
c  For usage:  in your par file
c     exactmodel = "multiBH"
c     kt_hubble       = 0.0   # 0.0 means MP solution
c     kt_nBH          = 2     # number of BHs (upto 4, currently)
c     m_bh1,co_bh1x,co_bh1y,co_bh1z =  1.0, -2.0,0.0,0.0  # masses and 
c     m_bh2,co_bh2x,co_bh2y,co_bh2z =  1.0,  2.0,0.0,0.0  #   locations
c-----------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine  multi_BHs(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx)

      implicit none

      DECLARE_CCTK_PARAMETERS

c     Input.
      CCTK_REAL x, y, z, t

c     Output.
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx

      logical firstcall

      CCTK_REAL kt_r, kt_aa, kt_omega
      CCTK_REAL one,zero
      CCTK_REAL kt_xbh(10),kt_ybh(10),kt_zbh(10),kt_mbh(10)

      integer i

      data firstcall /.true./
      save firstcall,kt_xbh,kt_ybh,kt_zbh,kt_mbh

c     Get parameters of the exact solution.

      if (firstcall) then

         write(*,*) ' welcome to Kastor-Traschen (Majumdar-Papapetrou)'

         if(kt_nBH.ge.1) then
            kt_xbh(1) = co_bh1x
            kt_ybh(1) = co_bh1y
            kt_zbh(1) = co_bh1z
            kt_mbh(1) = m_bh1
         endif

         if(kt_nBH.ge.2) then
            kt_xbh(2) = co_bh2x
            kt_ybh(2) = co_bh2y
            kt_zbh(2) = co_bh2z
            kt_mbh(2) = m_bh2
         endif

         if(kt_nBH.ge.3) then
            kt_xbh(3) = co_bh3x
            kt_ybh(3) = co_bh3y
            kt_zbh(3) = co_bh3z
            kt_mbh(3) = m_bh3
         endif

         if(kt_nBH.ge.4) then
            kt_xbh(4) = co_bh4x
            kt_ybh(4) = co_bh4y
            kt_zbh(4) = co_bh4z
            kt_mbh(4) = m_bh4
         endif

         write(*,*) 'time=',t
         write(*,*) '   mass    BH(x,y,z) '

         do i=1,kt_nBH
            write(*,'(4e12.3)') kt_mbh(i),kt_xbh(i),kt_ybh(i),kt_zbh(i)
         enddo

         firstcall = .false.

      end if

      one =1.0
      zero=0.0 

      kt_aa=exp(kt_hubble*t)
      kt_omega=1.0

      do i=1,kt_nBH
         kt_r=sqrt((x-kt_xbh(i))**2+(y-kt_ybh(i))**2+(z-kt_zbh(i))**2)
         kt_omega= kt_omega+kt_mbh(i)/kt_aa/kt_r
      enddo
      
c     write(*,*) kt_omega,kt_aa

      gdtt = -1.0/kt_omega**2
      gdtx = zero
      gdty = zero
      gdtz = zero
      gdxx = (kt_aa*kt_omega)**2
      gdyy = (kt_aa*kt_omega)**2
      gdzz = (kt_aa*kt_omega)**2
      gdxy = zero
      gdyz = zero
      gdzx = zero

      gutt = one/gdtt 
      gutx = zero
      guty = zero
      gutz = zero
      guxx = one/gdxx
      guyy = one/gdyy
      guzz = one/gdzz
      guxy = zero
      guyz = zero
      guzx = zero

      return
      end
