 /*@@
   @file      BrillLindquist.F
   @date      
   @author    
   @desc 
      Set up initial data for Brill Lindquist black holes
   @enddesc 
 @@*/

#include "cctk.h"
#include "declare_arguments.h"
#include "declare_parameters.h"

#include "../../../packages/CactusEinstein/Einstein/src/Einstein.h"

 /*@@
   @routine    BrillLindquist
   @date       
   @author     
   @desc 
      Set up initial data for Brill Lindquist black holes
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine BrillLindquist(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_PARAMETERS

      integer i
      REAL, dimension(4) :: hole_mass, hole_x0, hole_y0, hole_z0


c     Put parameters into arrays for following calculations
c     -----------------------------------------------------
      hole_x0(1)   = bl_x0_1
      hole_y0(1)   = bl_y0_1
      hole_z0(1)   = bl_z0_1
      hole_mass(1) = bl_M_1

      hole_x0(2)   = bl_x0_2
      hole_y0(2)   = bl_y0_2
      hole_z0(2)   = bl_z0_2
      hole_mass(2) = bl_M_2

      hole_x0(3)   = bl_x0_3
      hole_y0(3)   = bl_y0_3
      hole_z0(3)   = bl_z0_3
      hole_mass(3) = bl_M_3

      hole_x0(4)   = bl_x0_4
      hole_y0(4)   = bl_y0_4
      hole_z0(4)   = bl_z0_4
      hole_mass(4) = bl_M_4

c     Do not ask...
c     ----------
      hole_x0 = -hole_x0
      hole_y0 = -hole_y0
      hole_z0 = -hole_z0

c     Initialize to zero and then use +=
c     ----------------------------------
      psi = 1.0D0

      if (use_conformal_derivs == 1) then

        psix = 0.0D0
        psiy = 0.0D0
        psiz = 0.0D0
        psixx = 0.0D0
        psixy = 0.0D0
        psixz = 0.0D0
        psiyy = 0.0D0
        psiyz = 0.0D0
        psizz = 0.0D0

      end if

      do i = 1, bl_nbh

c     Maple Output
c     ------------
      Psi = Psi+hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+hole_x0(i)**2.D
     &0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*z*hole_z0(i)+h
     &ole_z0(i)**2.D0)/2.D0

      if (use_conformal_derivs == 1) then

      Psix = Psix-hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+hole_x0(i)**2
     &.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*z*hole_z0(i)
     &+hole_z0(i)**2.D0)**3.D0*(2.D0*x+2.D0*hole_x0(i))/4.D0
      Psiy = Psiy-hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+hole_x0(i)**2
     &.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*z*hole_z0(i)
     &+hole_z0(i)**2.D0)**3.D0*(2.D0*y+2.D0*hole_y0(i))/4.D0
      Psiz = Psiz-hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+hole_x0(i)**2
     &.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*z*hole_z0(i)
     &+hole_z0(i)**2.D0)**3.D0*(2.D0*z+2.D0*hole_z0(i))/4.D0
      Psixx = Psixx+3.D0/8.D0*hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+h
     &ole_x0(i)**2.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*
     &z*hole_z0(i)+hole_z0(i)**2.D0)**5.D0*(2.D0*x+2.D0*hole_x0(i))**2.D
     &0-hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+hole_x0(i)**2.D0+y**2+2
     &.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*z*hole_z0(i)+hole_z0(i
     &)**2.D0)**3.D0/2.D0
      Psixy = Psixy+3.D0/8.D0*hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+h
     &ole_x0(i)**2.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*
     &z*hole_z0(i)+hole_z0(i)**2.D0)**5.D0*(2.D0*x+2.D0*hole_x0(i))*(2.D
     &0*y+2.D0*hole_y0(i))
      Psixz = Psixz+3.D0/8.D0*hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+h
     &ole_x0(i)**2.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*
     &z*hole_z0(i)+hole_z0(i)**2.D0)**5.D0*(2.D0*x+2.D0*hole_x0(i))*(2.D
     &0*z+2.D0*hole_z0(i))
      Psiyy = Psiyy+3.D0/8.D0*hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+h
     &ole_x0(i)**2.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*
     &z*hole_z0(i)+hole_z0(i)**2.D0)**5.D0*(2.D0*y+2.D0*hole_y0(i))**2.D
     &0-hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+hole_x0(i)**2.D0+y**2+2
     &.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*z*hole_z0(i)+hole_z0(i
     &)**2.D0)**3.D0/2.D0
      Psiyz = Psiyz+3.D0/8.D0*hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+h
     &ole_x0(i)**2.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*
     &z*hole_z0(i)+hole_z0(i)**2.D0)**5.D0*(2.D0*y+2.D0*hole_y0(i))*(2.D
     &0*z+2.D0*hole_z0(i))
      Psizz = Psizz+3.D0/8.D0*hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+h
     &ole_x0(i)**2.D0+y**2+2.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*
     &z*hole_z0(i)+hole_z0(i)**2.D0)**5.D0*(2.D0*z+2.D0*hole_z0(i))**2.D
     &0-hole_mass(i)/sqrt(x**2+2.D0*x*hole_x0(i)+hole_x0(i)**2.D0+y**2+2
     &.D0*y*hole_y0(i)+hole_y0(i)**2.D0+z**2+2.D0*z*hole_z0(i)+hole_z0(i
     &)**2.D0)**3.D0/2.D0

      end if

      end do

c     Cactus conventions
c     ------------------
      if (use_conformal_derivs == 1) then

        psix = psix / psi
        psiy = psiy / psi
        psiz = psiz / psi
        psixx = psixx / psi
        psixy = psixy / psi
        psixz = psixz / psi
        psiyy = psiyy / psi
        psiyz = psiyz / psi
        psizz = psizz / psi

      end if

c     Metric depends on conformal choice
c     ----------------------------------
      if (conformal_state == CONFORMAL_METRIC) then

        gxx = 1d0
        gyy = 1d0
        gzz = 1d0
        gxy = 0d0 
        gxz = 0d0  
        gyz = 0d0   

      else 

        gxx = psi**4
        gyy = gxx
        gzz = gxx
        gxy = 0d0
        gyz = 0d0
        gzz = 0d0

      end if

c     Time-symmetric data
c     -------------------
      hxx = 0d0
      hyy = 0d0
      hzz = 0d0
      hxy = 0d0 
      hxz = 0d0  
      hyz = 0d0   

      return
      end

