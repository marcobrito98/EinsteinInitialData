 /*@@
   @file      Schwarzschild.c
   @date      Sun Oct 17 10:35:41 1999
   @author    Tom Goodale
   @desc
              C version of Scwhwarzschild lapse routine
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"

#include <string.h>

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Schwarzschild_c)

void Schwarzschild(CCTK_ARGUMENTS);


void Schwarzschild(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  const CCTK_REAL zero = 0.0, one = 1.0, two = 2.0, three = 3.0;
  CCTK_REAL tmp, r_squared, r_cubed;
  int i, npoints;

  npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  /*     conformal metric flag */
  if(CCTK_EQUALS(metric_type, "static conformal"))
  {
    int make_conformal_derivs;

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

    for (i = 0; i < npoints; i++)
    {
      /*        Compute conformal factor */
      psi[i] = ( one + mass/two/r[i]);

      if(make_conformal_derivs)
      {
        /*        derivatives of psi / psi */
        r_squared = r[i]*r[i];
        r_cubed   = r[i]*r_squared;
        tmp = mass/two/r_cubed/psi[i];


        psix[i] = -x[i]*tmp;
        psiy[i] = -y[i]*tmp;
        psiz[i] = -z[i]*tmp;

        if(*conformal_state > 2)
        {
          tmp = mass/two/(r_squared*r_cubed)/psi[i];
          psixy[i] = three*x[i]*y[i]*tmp;
          psixz[i] = three*x[i]*z[i]*tmp;
          psiyz[i] = three*y[i]*z[i]*tmp;
      
          psixx[i]  = (three*x[i]*x[i] - r_squared)*tmp;
          psiyy[i]  = (three*y[i]*y[i] - r_squared)*tmp;
          psizz[i]  = (three*z[i]*z[i] - r_squared)*tmp;
        }
      }
      gxx[i] = one;
      gyy[i] = one;
      gzz[i] = one;
      gxy[i] = zero;
      gxz[i] = zero;
      gyz[i] = zero;
    }
  }
  else
  {
    for (i = 0; i < npoints; i++)
    {
      r_squared = r[i] * r[i];
      gxx[i] = one + mass/two/(r_squared * r_squared);
      gyy[i] = gxx[i];
      gzz[i] = gxx[i];
      gxy[i] = zero;
      gxz[i] = zero;
      gyz[i] = zero;
    }
  }
  
  /*     If the initial lapse is not one ... */
  if (CCTK_Equals(initial_lapse,"schwarz"))
  {
    CCTK_INFO("Initialise with Schwarzschild lapse");

    for (i = 0; i < npoints; i++)
    {
      alp[i] = (2.*r[i] - mass)/(2.*r[i]+mass);
    }
  }

  /*     time symmetric initial slice */
  memset (kxx, 0, npoints * sizeof (CCTK_REAL));
  memset (kxy, 0, npoints * sizeof (CCTK_REAL));
  memset (kxz, 0, npoints * sizeof (CCTK_REAL));
  memset (kyy, 0, npoints * sizeof (CCTK_REAL));
  memset (kyz, 0, npoints * sizeof (CCTK_REAL));
  memset (kzz, 0, npoints * sizeof (CCTK_REAL));
}
