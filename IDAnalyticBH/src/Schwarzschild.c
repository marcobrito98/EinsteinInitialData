 /*@@
   @file      Schwarzschild.c
   @date      Sun Oct 17 10:35:41 1999
   @author    Tom Goodale
   @desc 
   C version of Scwhwarzschild lapse routine
   @enddesc 
 @@*/

#include <math.h>

#include "cctk.h"
#include "cctk_arguments.h"
#include "cctk_parameters.h"


/* Need include file from Einstein */
#include "CactusEinstein/Einstein/src/Einstein.h"

void Schwarzschild(CCTK_CARGUMENTS)
{
  DECLARE_CCTK_CARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL zero,one,two,three;
  CCTK_REAL tmp;
  CCTK_REAL r_squared;

  CCTK_REAL t1,r0;

  int i,j,k;
  int is,js,ks,ie,je,ke;
  CCTK_REAL inval;
  int nx,ny,nz;
  int index;

  zero = 0.0;
  one = 1.0;
  two = 2.0;
  three = 3.0;

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];
  
  /*     conformal metric flag */
  if(use_conformal == 1)
  {
  
    *conformal_state = CONFORMAL_METRIC;


    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);

          /*        Compute conformal factor */
          psi[index] = ( one + mass/two/r[index]);
         

          /*        derivatives of psi / psi */

          tmp = mass/two/pow(r[index],3) / psi[index];
          psix[index] = -x[index]*tmp;
          psiy[index] = -y[index]*tmp;
          psiz[index] = -z[index]*tmp;
         
          tmp = mass/two/pow(r[index],5)/psi[index];
          psixy[index] = three*x[index]*y[index]*tmp;
          psixz[index] = three*x[index]*z[index]*tmp;
          psiyz[index] = three*y[index]*z[index]*tmp;

          r_squared = r[index]*r[index];
          psixx[index]  = (three*x[index]*x[index] - r_squared)*tmp;
          psiyy[index]  = (three*y[index]*y[index] - r_squared)*tmp;
          psizz[index]  = (three*z[index]*z[index] - r_squared)*tmp;
         
          gxx[index] = one;
          gyy[index] = one;
          gzz[index] = one;
          gxy[index] = zero; 
          gxz[index] = zero;  
          gyz[index] = zero;  
        }
      }
    }
  }     
  else
  {
    *conformal_state = NOCONFORMAL_METRIC;

    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);

          gxx[index] = pow(( one + mass/two/r[index]), 4);
          gyy[index] = gxx[index];
          gzz[index] = gxx[index];
          gxy[index] = zero;
          gxz[index] = zero;
          gyz[index] = zero;
        }
      }
    }
  }

  /*     If the initial lapse is not one ... */
  if (CCTK_Equals(initial_lapse,"schwarz"))
  { 
    CCTK_INFO("Initialise with schwarzschild lapse");

    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);

          alp[index] = (2.*r[index] - mass)/(2.*r[index]+mass);
        }
      }
    }
  }
      
  /*     time symmetric initial slice */

  for(k=0; k < nz; k++)
  {
    for(j=0; j < ny; j++)
    {
      for(i=0; i < nx; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH, i,j,k);
        
        kxx[index] = zero;
        kxy[index] = zero;
        kxz[index] = zero;
        kyy[index] = zero;
        kyz[index] = zero;
        kzz[index] = zero;
      }
    }
  }
    
  return;
}
