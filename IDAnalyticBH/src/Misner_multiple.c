 /*@@
   @file      Misner_multiple.F
   @date      
   @author    Carsten Gundlach
   @desc 
              Set up initial data for multiple Misner black holes
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"

#include <string.h>

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Misner_multiple_c)

#define SQR(a) ((a)*(a))

void Misner_init(int n, CCTK_REAL mu, int terms);
void MisnerEvalPsi(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *res); 
void Misner_multiple(CCTK_ARGUMENTS);

 /*@@
   @routine    Misner_multiple
   @date       
   @author     Carsten Gundlach
   @desc 
               Set up initial data for multiple Misner black holes
   @enddesc 
   @calls      MisnerEvalPsi
   @hdate Fri Apr 26 10:04:05 2002 @hauthor Tom Goodale
   @hdesc Changed to use new StaticConformal stuff
   @endhistory 
@@*/
void Misner_multiple(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i, npoints;
  CCTK_REAL xval, yval, zval;
  CCTK_REAL inv_psi, tmp0, tmp1, tmp2, tmp3, tmp4;
  const CCTK_REAL nm_eps             = 1e-6,   /* finite differencing step*/
                  halved_inv_nm_eps  = 5e+5,
                  inv_nm_eps_squared = 1e+12;
  const CCTK_REAL one = 1.0;
  int make_conformal_derivs;


  /* Check if we should create and store conformal factor stuff */
  if(CCTK_EQUALS(metric_type, "static conformal"))
  {
    *conformal_state = 1;

    if(CCTK_EQUALS(conformal_storage,"factor+derivs"))
    {
      *conformal_state = 2;
      make_conformal_derivs = 1;
    }
    else if(CCTK_EQUALS(conformal_storage,"factor+derivs+2nd derivs"))
    {
      *conformal_state = 3;
      make_conformal_derivs = 1;
    }
  }      
  else
  {
    make_conformal_derivs = 0;
  }


  npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  /*     Initialize C global variables
   *     -----------------------------
   */

  Misner_init(misner_nbh, mu, nmax);

  /*     Get value of psi pointwise from a C function
   *     --------------------------------------------
   */

  for(i = 0; i < npoints; i++)
  {
    xval = x[i];
    yval = y[i];
    zval = z[i];

    MisnerEvalPsi(xval, yval, zval, &tmp0);
    psi[i] = tmp0;
           
    /*              Only calculate derivatives of psi if required
     *              ---------------------------------------------
     */
    if (make_conformal_derivs)
    {
      MisnerEvalPsi(xval+nm_eps,yval,zval,&tmp1);
      MisnerEvalPsi(xval-nm_eps,yval,zval,&tmp2);
      psix[i] = (tmp1-tmp2) * halved_inv_nm_eps;

      if(*conformal_state > 2)
      {
        psixx[i] = (tmp1+tmp2-2.0*tmp0) * inv_nm_eps_squared;
      }

      MisnerEvalPsi(xval,yval+nm_eps,zval,&tmp1);
      MisnerEvalPsi(xval,yval-nm_eps,zval,&tmp2);
      psiy[i]  = (tmp1-tmp2) * halved_inv_nm_eps;

      if(*conformal_state > 2)
      {
        psiyy[i] = (tmp1+tmp2-2.0*tmp0) * inv_nm_eps_squared;
      }

      MisnerEvalPsi(xval,yval,zval+nm_eps,&tmp1);
      MisnerEvalPsi(xval,yval,zval-nm_eps,&tmp2);
      psiz[i] = (tmp1-tmp2) * halved_inv_nm_eps;

      if(*conformal_state > 2)
      {
        psizz[i] = (tmp1+tmp2-2.0*tmp0) * inv_nm_eps_squared;
      }

      if(*conformal_state > 2)
      {
        MisnerEvalPsi(xval+nm_eps,yval+nm_eps,zval,&tmp1);
        MisnerEvalPsi(xval+nm_eps,yval-nm_eps,zval,&tmp2);
        MisnerEvalPsi(xval-nm_eps,yval+nm_eps,zval,&tmp3);
        MisnerEvalPsi(xval-nm_eps,yval-nm_eps,zval,&tmp4);
        psixy[i] = 0.25*(tmp1-tmp2-tmp3+tmp4) * inv_nm_eps_squared;

        MisnerEvalPsi(xval,yval+nm_eps,zval+nm_eps,&tmp1);
        MisnerEvalPsi(xval,yval-nm_eps,zval+nm_eps,&tmp2);
        MisnerEvalPsi(xval,yval+nm_eps,zval-nm_eps,&tmp3);
        MisnerEvalPsi(xval,yval-nm_eps,zval-nm_eps,&tmp4);
        psiyz[i] = 0.25*(tmp1-tmp2-tmp3+tmp4) * inv_nm_eps_squared;

        MisnerEvalPsi(xval+nm_eps,yval,zval+nm_eps,&tmp1);
        MisnerEvalPsi(xval+nm_eps,yval,zval-nm_eps,&tmp2);
        MisnerEvalPsi(xval-nm_eps,yval,zval+nm_eps,&tmp3);
        MisnerEvalPsi(xval-nm_eps,yval,zval-nm_eps,&tmp4);
        psixz[i] = 0.25*(tmp1-tmp2-tmp3+tmp4) * inv_nm_eps_squared;
      }
    }
  }

  /*     Cactus conventions
   *     ------------------
   */

  if (make_conformal_derivs)
  {
    for(i = 0; i < npoints; i++)
    {
      inv_psi = one / psi[i];

      psix[i]  *= inv_psi;
      psiy[i]  *= inv_psi;
      psiz[i]  *= inv_psi;
      if(*conformal_state > 2)
      {
        psixx[i] *= inv_psi;
        psixy[i] *= inv_psi;
        psixz[i] *= inv_psi;
        psiyy[i] *= inv_psi;
        psiyz[i] *= inv_psi;
        psizz[i] *= inv_psi;
      }
    }
  }

  /*     Metric depends on conformal state
   *     ---------------------------------
   */

  if (CCTK_EQUALS(metric_type, "static conformal"))
  {
    for(i = 0; i < npoints; i++)
    {
      gxx[i] = one;
      gyy[i] = one;
      gzz[i] = one;
    }
  }
  else 
  {
    for(i = 0; i < npoints; i++)
    {
      gxx[i] = psi[i] * psi[i] * psi[i] * psi[i];
      gyy[i] = gxx[i];
      gzz[i] = gxx[i];
    }
  }
  memset (gxy, 0, npoints * sizeof (gxy [0]));
  memset (gxz, 0, npoints * sizeof (gxz [0]));
  memset (gyz, 0, npoints * sizeof (gyz [0]));

  /*     Time-symmetric data
   *     -------------------
   */

  memset (kxx, 0, npoints * sizeof (kxx [0]));
  memset (kyy, 0, npoints * sizeof (kyy [0]));
  memset (kzz, 0, npoints * sizeof (kzz [0]));
  memset (kxy, 0, npoints * sizeof (kxy [0]));
  memset (kxz, 0, npoints * sizeof (kxz [0]));
  memset (kyz, 0, npoints * sizeof (kyz [0]));

}
