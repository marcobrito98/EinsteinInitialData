 /*@@
   @file      Misner_standard.c
   @date      March 1997
   @author    Joan Masso
   @desc
              Set up initial data for two Misner black holes
   @enddesc
   @history
   @hdate Sun Oct 17 11:05:48 1999 @hauthor Tom Goodale
   @hdesc Converted to C
   @endhistory
   @version   $Id$
 @@*/

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Misner_standard_c)

#define SQR(a) ((a)*(a))

void Misner_standard(CCTK_ARGUMENTS);

 /*@@
   @routine    Misner_standard
   @date
   @author     Joan Masso, Ed Seidel
   @desc
            Initialize the metric with a time symmetrical
            black hole spacetime containing
            two axially symmetric misner black holes with a
            mass/length parameter mu. The mass is computed.
            The spacetime line element has the form:
                 $$ ds^2 = -dt^2 + \Psi^4 (dx^2+dy^2+dz^2) $$
            and only $\Psi$ differs.
            (Conformal factor from Karen Camarda)
   @enddesc

   @par      mu
   @pdesc    Misner parameter.
   @ptype    real
   @pcomment Values less than 1.8 do not really correspond to two
            black holes, as there is an initial single event horizon
            surrounding the throats. So, with low values of mu we also
            have distorted single black holes.
   @endpar

   @par      nmax
   @pdesc    Summation limit for the misner series in the 'twobh' case.
   @ptype    integer
   @endpar

@@*/
void Misner_standard(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i, n;
  int npoints;
  CCTK_REAL *csch, *coth, inv_r1, inv_r2;
  CCTK_REAL x_squared, y_squared, xy_squared;
  CCTK_REAL inv_r1_cubed, inv_r2_cubed;
  CCTK_REAL inv_r1_5, inv_r2_5;
  CCTK_REAL inv_psi;
  CCTK_INT powfac;
  CCTK_INT adm_mass;
  const CCTK_REAL zero = 0.0, one = 1.0, three = 3.0;


  /* total number of points on this processor */
  npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];


  /*     Initialize so we can accumulate
   *     -------------------------------
   */
  if (use_conformal_derivs)
  {
    memset (psix, 0, npoints * sizeof (CCTK_REAL));
    memset (psiy, 0, npoints * sizeof (CCTK_REAL));
    memset (psiz, 0, npoints * sizeof (CCTK_REAL));
    memset (psixx, 0, npoints * sizeof (CCTK_REAL));
    memset (psixy, 0, npoints * sizeof (CCTK_REAL));
    memset (psixz, 0, npoints * sizeof (CCTK_REAL));
    memset (psiyy, 0, npoints * sizeof (CCTK_REAL));
    memset (psiyz, 0, npoints * sizeof (CCTK_REAL));
    memset (psizz, 0, npoints * sizeof (CCTK_REAL));
  }

  csch = (CCTK_REAL *) malloc (2 * (nmax + 1) * sizeof (CCTK_REAL));
  coth = csch + nmax + 1;

  /*     compute the ADM mass
   *     --------------------
   */
  adm_mass = zero;
  for(n = 1; n <= nmax; n++)
  {
    csch[n] = one / sinh(mu*n);
    coth[n] = one / tanh(mu*n);
    adm_mass   += 4.0 * csch[n];
  }
  CCTK_VInfo(CCTK_THORNSTRING, "ADM mass is %f", (double) adm_mass);


  for(i = 0; i < npoints; i++)
  {
    psi [i] = one;

    x_squared  = SQR(x[i]);
    y_squared  = SQR(y[i]);
    xy_squared = x_squared + y_squared;

    for(n = nmax; n >= 1; n--)
    {
      inv_r1 = one / sqrt(xy_squared+SQR(z[i]+coth[n]));
      inv_r2 = one / sqrt(xy_squared+SQR(z[i]-coth[n]));

      psi[i] += csch[n]*(inv_r1 + inv_r2);

      if (use_conformal_derivs)
      {
        inv_r1_cubed = inv_r1 * inv_r1 * inv_r1;
        inv_r2_cubed = inv_r2 * inv_r2 * inv_r2;
        inv_r1_5     = inv_r1 * inv_r1 * inv_r1_cubed;
        inv_r2_5     = inv_r2 * inv_r2 * inv_r2_cubed;
        psix[i]  +=  -x[i] * (inv_r2_cubed + inv_r1_cubed) * csch[n];
        psiy[i]  +=  -y[i] * (inv_r2_cubed + inv_r1_cubed) * csch[n];
        psiz[i]  +=  (-(z[i]-coth[n])*inv_r2_cubed - (z[i]+coth[n])*inv_r1_cubed) * csch[n];
        psixx[i] +=  (three*x_squared*(inv_r1_5 + inv_r2_5)
                      - inv_r1_cubed - inv_r2_cubed) * csch[n];
        psixy[i] +=  three*x[i]*y[i]*(inv_r1_5 + inv_r2_5) * csch[n];
        psixz[i] +=  (three*x[i]*(z[i] - coth[n])*inv_r2_5
                      + three*x[i]*(z[i] + coth[n])*inv_r1_5) * csch[n];
        psiyy[i] +=  (three*y_squared*(inv_r1_5 + inv_r2_5)
                      - inv_r1_cubed - inv_r2_cubed) * csch[n];
        psiyz[i] +=  (three*y[i]*(z[i] - coth[n])*inv_r2_5
                      + three*y[i]*(z[i] + coth[n])*inv_r1_5) * csch[n];
        psizz[i] += (-inv_r2_cubed+three*SQR(z[i] - coth[n])*inv_r2_5
                      + three*SQR(z[i] + coth[n])*inv_r1_5 - inv_r1_cubed) * csch[n];
      }
    }

    /*     Cactus convention
     *     -----------------
     */
    if (use_conformal_derivs)
    {
      inv_psi = one / psi[i];

      psix[i]  *= inv_psi;
      psiy[i]  *= inv_psi;
      psiz[i]  *= inv_psi;
      psixx[i] *= inv_psi;
      psixy[i] *= inv_psi;
      psixz[i] *= inv_psi;
      psiyy[i] *= inv_psi;
      psiyz[i] *= inv_psi;
      psizz[i] *= inv_psi;
    }
  }

  /*     Should initialize lapse to Cadez value if possible
   *     --------------------------------------------------
   */

  if (CCTK_Equals(initial_lapse,"cadez"))
  {
    CCTK_INFO("Initialise with cadez lapse");

    for(i = 0; i < npoints; i++)
    {
      xy_squared = SQR(x[i]) + SQR(y[i]);

      alp[i] = one;

      powfac = 1;

      for(n = 1; n <= nmax; n++)
      {
        inv_r1 = one / sqrt(xy_squared+SQR(z[i]+coth[n]));
        inv_r2 = one / sqrt(xy_squared+SQR(z[i]-coth[n]));
        powfac = -powfac;

        alp[i] += powfac * csch[n] * (inv_r1 + inv_r2);
      }

      alp[i] /= psi[i];
    }
  }

  /*     Metric depends on conformal state
   *     ---------------------------------
   */

  if (*conformal_state == CONFORMAL_METRIC)
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
  memset (gxy, 0, npoints * sizeof (CCTK_REAL));
  memset (gxz, 0, npoints * sizeof (CCTK_REAL));
  memset (gyz, 0, npoints * sizeof (CCTK_REAL));

  /*     Time-symmetric data
   *     -------------------
   */
  memset (kxx, 0, npoints * sizeof (CCTK_REAL));
  memset (kyy, 0, npoints * sizeof (CCTK_REAL));
  memset (kzz, 0, npoints * sizeof (CCTK_REAL));
  memset (kxy, 0, npoints * sizeof (CCTK_REAL));
  memset (kxz, 0, npoints * sizeof (CCTK_REAL));
  memset (kyz, 0, npoints * sizeof (CCTK_REAL));

  if (csch)
  {
    free (csch);
  }
}
