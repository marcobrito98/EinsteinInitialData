 /*@@
   @file      BrillLindquist.F
   @date
   @author
   @desc
              Set up initial data for Brill Lindquist black holes
   @enddesc
   @version   $Header$
 @@*/

#include "cctk.h"

#include <math.h>
#include <string.h>

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_BrillLindquist_c)


#define SQR(a) ((a)*(a))

#define MAX_HOLES 4

void BrillLindquist(CCTK_ARGUMENTS);

 /*@@
   @routine    BrillLindquist
   @date
   @author
   @desc
               Set up initial data for Brill Lindquist black holes
   @enddesc
   @hdate Fri Apr 26 10:04:05 2002 @hauthor Tom Goodale
   @hdesc Changed to use new StaticConformal stuff
   @endhistory 
@@*/
void BrillLindquist(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int n;
  CCTK_REAL hole_mass[MAX_HOLES], hole_x0[MAX_HOLES];
  CCTK_REAL hole_y0[MAX_HOLES], hole_z0[MAX_HOLES];
  CCTK_REAL tmp1, tmp2, tmp3;
  CCTK_REAL xval, yval, zval;
  CCTK_REAL x_2, y_2, z_2;
  int i, npoints;
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

  
  if (make_conformal_derivs == 1)
  {
    memset (psix, 0, npoints * sizeof (CCTK_REAL));
    memset (psiy, 0, npoints * sizeof (CCTK_REAL));
    memset (psiz, 0, npoints * sizeof (CCTK_REAL));

    if(*conformal_state > 2)
    {
      memset (psixx, 0, npoints * sizeof (CCTK_REAL));
      memset (psixy, 0, npoints * sizeof (CCTK_REAL));
      memset (psixz, 0, npoints * sizeof (CCTK_REAL));
      memset (psiyy, 0, npoints * sizeof (CCTK_REAL));
      memset (psiyz, 0, npoints * sizeof (CCTK_REAL));
      memset (psizz, 0, npoints * sizeof (CCTK_REAL));
    }
  }

  for (i = 0; i < npoints; i++)
  {
    /*     Initialize to zero and then use +=
     *     ----------------------------------
     */

    psi[i] = 1.0;

    x_2 = SQR(x[i]);
    y_2 = SQR(y[i]);
    z_2 = SQR(z[i]);

    xval = x[i];
    yval = y[i];
    zval = z[i];

    for(n = 0;  n < bl_nbh; n++)
    {

      /*     Maple Output
       *     ------------
       */

      tmp1 = sqrt(x_2+2.0*xval*hole_x0[n]+SQR(hole_x0[n])
                 +y_2+2.0*yval*hole_y0[n]+SQR(hole_y0[n])
                 +z_2+2.0*zval*hole_z0[n]+SQR(hole_z0[n]));

      psi[i] += hole_mass[n]/tmp1*0.5;

      if (make_conformal_derivs == 1)
      {
        tmp2 = 1 / (tmp1 * tmp1 * tmp1);
        tmp3 = 4 * (3.0 / 8.0) * hole_mass[n] * tmp2 / (tmp1 * tmp1);
        tmp2 *= -0.5 * hole_mass[n];

        psix[i] +=  tmp2 * (xval+hole_x0[n]);
        psiy[i] +=  tmp2 * (yval+hole_y0[n]);
        psiz[i] +=  tmp2 * (zval+hole_z0[n]);

        if(*conformal_state > 2)
        {
          psixx[i] += tmp3 * SQR(xval+hole_x0[n]) + tmp2;
          
          psixy[i] += tmp3 * (xval+hole_x0[n]) * (yval+hole_y0[n]);
          psixz[i] += tmp3 * (xval+hole_x0[n]) * (zval+hole_z0[n]);
          
          psiyy[i] += tmp3 * SQR(yval+hole_y0[n]) + tmp2;
          psiyz[i] += tmp3 * (yval+hole_y0[n]) * (zval+hole_z0[n]);
          
          psizz[i] += tmp3 * SQR(zval+hole_z0[n]) + tmp2;
        }
      }
    }
  }

  /*     Cactus conventions
   *     ------------------
   */

  if (make_conformal_derivs == 1)
  {
    for (i = 0; i < npoints; i++)
    {
      tmp1 = 1 / psi[i];

      psix[i]  *= tmp1;
      psiy[i]  *= tmp1;
      psiz[i]  *= tmp1;

      if(*conformal_state > 2)
      {
        psixx[i] *= tmp1;
        psixy[i] *= tmp1;
        psixz[i] *= tmp1;
        psiyy[i] *= tmp1;
        psiyz[i] *= tmp1;
        psizz[i] *= tmp1;
      }
    }
  }

  /*     Metric depends on conformal state
   *     ---------------------------------
   */

  if (CCTK_EQUALS(metric_type, "static conformal"))
  {
    for (i = 0; i < npoints; i++)
    {
      gxx[i] = 1.0;
      gyy[i] = 1.0;
      gzz[i] = 1.0;
      gxy[i] = 0.0;
      gxz[i] = 0.0;
      gyz[i] = 0.0;
    }
  }
  else
  {
    for (i = 0; i < npoints; i++)
    {
      tmp1 = psi[i] * psi[i];
      gxx[i] = tmp1 * tmp1;
      gyy[i] = gxx[i];
      gzz[i] = gxx[i];
      gxy[i] = 0.0;
      gxz[i] = 0.0;
      gyz[i] = 0.0;
    }
  }

  /*     Time-symmetric data
   *     -------------------
   */
  memset (kxx, 0, npoints * sizeof (CCTK_REAL));
  memset (kyy, 0, npoints * sizeof (CCTK_REAL));
  memset (kzz, 0, npoints * sizeof (CCTK_REAL));
  memset (kxy, 0, npoints * sizeof (CCTK_REAL));
  memset (kxz, 0, npoints * sizeof (CCTK_REAL));
  memset (kyz, 0, npoints * sizeof (CCTK_REAL));
}
