 /*@@
   @file      Misner_multiple.F
   @date      
   @author    Carsten Gundlach
   @desc 
      Set up initial data for multiple Misner black holes
   @enddesc 
 @@*/

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

#define SQR(a) ((a)*(a))

void Misner_init(int n, CCTK_REAL mu, int terms);
void MisnerEvalPsi(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *res); 

 /*@@
   @routine    Misner_multiple
   @date       
   @author     Carsten Gundlach
   @desc 
      Set up initial data for multiple Misner black holes
   @enddesc 
   @calls     
   @history 
 
   @endhistory 

@@*/
void Misner_multiple(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i, j, k;
  CCTK_REAL xval, yval, zval;
  CCTK_REAL one, zero, tmp0, tmp1, tmp2, tmp3, tmp4;
  CCTK_REAL nm_eps = 0.000001;   /* finite differencing step*/

  int nx,ny,nz;
  int index;

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];

  one  = 1.0;
  zero = 0.0;

  /*     Initialize C global variables
   *     -----------------------------
   */

  Misner_init(misner_nbh,mu,nmax);

  /*     Get value of psi pointwise from a C function
   *     --------------------------------------------
   */
  for(k=0; k < nz; k++)
  {
    for(j=0; j < ny; j++)
    {
      for(i=0; i < nx; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH, i,j,k);

        xval = x[index];
        yval = y[index];
        zval = z[index];

        MisnerEvalPsi(xval,yval,zval,&tmp0);
        psi[index] = tmp0;
               
        /*              Only calculate derivatives of psi if required
         *              ---------------------------------------------
         */
        if (use_conformal_derivs == 1)
        {
          MisnerEvalPsi(xval+nm_eps,yval,zval,&tmp1);
          MisnerEvalPsi(xval-nm_eps,yval,zval,&tmp2);
          psix[index] = 0.5*(tmp1-tmp2)/nm_eps;
          psixx[index] = (tmp1+tmp2-2.0*tmp0)/SQR(nm_eps);

          MisnerEvalPsi(xval,yval+nm_eps,zval,&tmp1);
          MisnerEvalPsi(xval,yval-nm_eps,zval,&tmp2);
          psiy[index]  = 0.5*(tmp1-tmp2)/nm_eps;
          psiyy[index] = (tmp1+tmp2-2.0*tmp0)/SQR(nm_eps);

          MisnerEvalPsi(xval,yval,zval+nm_eps,&tmp1);
          MisnerEvalPsi(xval,yval,zval-nm_eps,&tmp2);
          psiz[index] = 0.5*(tmp1-tmp2)/nm_eps;
          psizz[index] = (tmp1+tmp2-2.0*tmp0)/SQR(nm_eps);

          MisnerEvalPsi(xval+nm_eps,yval+nm_eps,zval,&tmp1);
          MisnerEvalPsi(xval+nm_eps,yval-nm_eps,zval,&tmp2);
          MisnerEvalPsi(xval-nm_eps,yval+nm_eps,zval,&tmp3);
          MisnerEvalPsi(xval-nm_eps,yval-nm_eps,zval,&tmp4);
          psixy[index] = 0.25*(tmp1-tmp2-tmp3+tmp4)/SQR(nm_eps);

          MisnerEvalPsi(xval,yval+nm_eps,zval+nm_eps,&tmp1);
          MisnerEvalPsi(xval,yval-nm_eps,zval+nm_eps,&tmp2);
          MisnerEvalPsi(xval,yval+nm_eps,zval-nm_eps,&tmp3);
          MisnerEvalPsi(xval,yval-nm_eps,zval-nm_eps,&tmp4);
          psiyz[index] = 0.25*(tmp1-tmp2-tmp3+tmp4)/SQR(nm_eps);

          MisnerEvalPsi(xval+nm_eps,yval,zval+nm_eps,&tmp1);
          MisnerEvalPsi(xval+nm_eps,yval,zval-nm_eps,&tmp2);
          MisnerEvalPsi(xval-nm_eps,yval,zval+nm_eps,&tmp3);
          MisnerEvalPsi(xval-nm_eps,yval,zval-nm_eps,&tmp4);
          psixz[index] = 0.25*(tmp1-tmp2-tmp3+tmp4)/SQR(nm_eps);

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
          gxy[index] = zero;
          gxz[index] = zero;
          gyz[index] = zero;
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
          
        kxx[index] = zero;
        kyy[index] = zero;
        kzz[index] = zero;
        kxy[index] = zero;
        kxz[index] = zero;
        kyz[index] = zero;
      }
    }
  }

  return;  

}
