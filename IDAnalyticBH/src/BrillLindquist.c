 /*@@
   @file      BrillLindquist.F
   @date      
   @author    
   @desc 
      Set up initial data for Brill Lindquist black holes
   @enddesc 
 @@*/

#include <math.h>

#include "cctk.h"
#include "cctk_arguments.h"
#include "cctk_parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

#define SQR(a) ((a)*(a))

#define MAX_HOLES 4
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

void BrillLindquist(CCTK_CARGUMENTS)
{
  DECLARE_CCTK_CARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int n;
  CCTK_REAL hole_mass[MAX_HOLES], hole_x0[MAX_HOLES], hole_y0[MAX_HOLES], hole_z0[MAX_HOLES];
  CCTK_REAL tmp1, tmp2, tmp3;
  CCTK_REAL xval, yval, zval;
  CCTK_REAL x_2, y_2, z_2;
  int i,j,k;
  int nx,ny,nz;
  int index;

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];

  /*     Put parameters into arrays for following calculations
   *     -----------------------------------------------------
   */
  hole_x0[0]   = -bl_x0_1;
  hole_y0[0]   = -bl_y0_1;
  hole_z0[0]   = -bl_z0_1;
  hole_mass[0] = bl_M_1;

  hole_x0[1]   = -bl_x0_2;
  hole_y0[1]   = -bl_y0_2;
  hole_z0[1]   = -bl_z0_2;
  hole_mass[1] = bl_M_2;

  hole_x0[2]   = -bl_x0_3;
  hole_y0[2]   = -bl_y0_3;
  hole_z0[2]   = -bl_z0_3;
  hole_mass[2] = bl_M_3;

  hole_x0[3]   = -bl_x0_4;
  hole_y0[3]   = -bl_y0_4;
  hole_z0[3]   = -bl_z0_4;
  hole_mass[3] = bl_M_4;


  for(k=0; k < nz; k++)
  {
    for(j=0; j < ny; j++)
    {
      for(i=0; i < nx; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          

        /*     Initialize to zero and then use +=
         *     ----------------------------------
         */

        psi[index] = 1.0;

        if (use_conformal_derivs == 1)
        {
          psix[index]  = 0.0;
          psiy[index]  = 0.0;
          psiz[index]  = 0.0;
          psixx[index] = 0.0;
          psixy[index] = 0.0;
          psixz[index] = 0.0;
          psiyy[index] = 0.0;
          psiyz[index] = 0.0;
          psizz[index] = 0.0;
        }

        x_2 = SQR(x[index]);
        y_2 = SQR(y[index]);

        xval = x[index];
        yval = y[index];
        zval = z[index];

        for(n = 0;  n < bl_nbh; i++)
        {

          /*     Maple Output
           *     ------------
           */
          
          tmp1 = sqrt(x_2+2.0*xval*hole_x0[n]+SQR(hole_x0[n])
                     +y_2+2.0*yval*hole_y0[n]+SQR(hole_y0[n])
                     +z_2+2.0*zval*hole_z0[n]+SQR(hole_z0[n]));

          psi[index] += hole_mass[n]/tmp1/2.0;

          if (use_conformal_derivs == 1)
          {
            tmp2 = pow(tmp1, 3);

            psix[index] +=  -hole_mass[n]/tmp2*(2.0*xval+2.0*hole_x0[n])/4.0;
            psiy[index] +=  -hole_mass[n]/tmp2*(2.0*yval+2.0*hole_y0[n])/4.0;
            psiz[index] +=  -hole_mass[n]/tmp2*(2.0*zval+2.0*hole_z0[n])/4.0;

            tmp2 = pow(tmp1, 5.0);
            tmp3 = pow(tmp1, 3.0);

            psixx[index] += 3.0/8.0*hole_mass[n]/tmp2*SQR(2.0*xval+2.0*hole_x0[n])-hole_mass[n]/tmp3/2.0;
        
            psixy[index] += 3.0/8.0*hole_mass[n]/tmp2*(2.0*xval+2.0*hole_x0[n])*(2.0*yval+2.0*hole_y0[n]);
            psixz[index] += 3.0/8.0*hole_mass[n]/tmp2*(2.0*xval+2.0*hole_x0[n])*(2.0*zval+2.0*hole_z0[n]);

            psiyy[index] += 3.0/8.0*hole_mass[n]/tmp2*SQR(2.0*yval+2.0*hole_y0[n])-hole_mass[n]/tmp3/2.0;
            psiyz[index] += 3.0/8.0*hole_mass[n]/tmp2*(2.0*yval+2.0*hole_y0[n])*(2.0*zval+2.0*hole_z0[n]);

            psizz[index] += 3.0/8.0*hole_mass[n]/tmp2*SQR(2.0*zval+2.0*hole_z0[n])-hole_mass[n]/tmp3/2.0;
          }
        }
      }
    }
  }

  /*     Cactus conventions
   *     ------------------
   */

  if (use_conformal_derivs == 1)
  {
    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          
          psix[index]  /= psi[index];
          psiy[index]  /= psi[index];
          psiz[index]  /= psi[index];
          psixx[index] /= psi[index];
          psixy[index] /= psi[index];
          psixz[index] /= psi[index];
          psiyy[index] /= psi[index];
          psiyz[index] /= psi[index];
          psizz[index] /= psi[index];
        }
      }
    }
  }

  /*     Metric depends on conformal state
   *     ---------------------------------
   */

  if (*conformal_state == CONFORMAL_METRIC)
  {
    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          

          gxx[index] = 1.0;
          gyy[index] = 1.0;
          gzz[index] = 1.0;
          gxy[index] = 0.0;
          gxz[index] = 0.0;
          gyz[index] = 0.0;
        }
      }
    }
  }
  else 
  {
    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          
          gxx[index] = pow(psi[index],4);
          gyy[index] = gxx[index];
          gzz[index] = gxx[index];
          gxy[index] = 0.0;
          gxz[index] = 0.0;
          gyz[index] = 0.0;
        }
      }
    }
  }

  /*     Time-symmetric data
   *     -------------------
   */

  for(k=0; k < nz; k++)
  {
    for(j=0; j < ny; j++)
    {
      for(i=0; i < nx; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          
        kxx[index] = 0.0;
        kyy[index] = 0.0;
        kzz[index] = 0.0;
        kxy[index] = 0.0;
        kxz[index] = 0.0;
        kyz[index] = 0.0;
      }
    }
  }

  return;
}
