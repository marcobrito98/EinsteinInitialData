// TwoPunctures:  File  "TwoPunctures.c"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

static inline double pow2 (const double x)
{
  return x*x;
}

static inline double pow4 (const double x)
{
  return x*x*x*x;
}

// -------------------------------------------------------------------
void
TwoPunctures (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int i, j, k, ntotal = n1 * n2 * n3 * nvar;
  double *F;
  derivs u, v;

  F = dvector (0, ntotal - 1);
  allocate_derivs (&u, ntotal);
  allocate_derivs (&v, ntotal);

  CCTK_INFO ("Solving puncture equation");
  Newton (nvar, n1, n2, n3, v, 1.e-10, 5);

  F_of_v (nvar, n1, n2, n3, v, F, u);

  CCTK_INFO ("Interpolating result");
  if (CCTK_EQUALS(metric_type, "static conformal")) {
    if (CCTK_EQUALS(conformal_storage, "factor")) {
      *conformal_state = 1;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs")) {
      *conformal_state = 2;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs")) {
      *conformal_state = 3;
    }
  } else {
    *conformal_state = 0;
  }

  for (k = 0; k < cctk_lsh[2]; ++k)
  {
    for (j = 0; j < cctk_lsh[1]; ++j)
    {
      for (i = 0; i < cctk_lsh[0]; ++i)
      {

	const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
        
        const double r_plus
          = sqrt(pow2(x[ind] - par_b) + pow2(y[ind]) + pow2(z[ind]));
        const double r_minus
          = sqrt(pow2(x[ind] + par_b) + pow2(y[ind]) + pow2(z[ind]));
        
        const double U = PunctIntPolAtArbitPosition
	  (0, nvar, n1, n2, n3, v, x[ind], y[ind], z[ind]);
        const double psi1 = 1
          + 0.5 * par_m_plus / r_plus
          + 0.5 * par_m_minus / r_minus + U;
        double static_psi = 1;
        
        double Aij[3][3];
        BY_Aijofxyz (x[ind], y[ind], z[ind], Aij);

        if (*conformal_state > 0) {

          double xp, yp, zp, rp, ir;
          double s1, s3, s5;
          double p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;

          p = 1.0;
          px = py = pz = 0.0;
          pxx = pxy = pxz = 0.0;
          pyy = pyz = pzz = 0.0;

          /* first puncture */
          xp = x[ind] - par_b;
          yp = y[ind];
          zp = z[ind];
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          ir = 1.0/rp;

          s1 = 0.5*par_m_plus*ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          /* second puncture */
          xp = x[ind] + par_b;
          yp = y[ind];
          zp = z[ind];
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          ir = 1.0/rp;

          s1 = 0.5*par_m_minus*ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          if (*conformal_state >= 1) {
            static_psi = p;
            psi[ind] = static_psi;
          }
          if (*conformal_state >= 2) {
            psix[ind] = px / static_psi;
            psiy[ind] = py / static_psi;
            psiz[ind] = pz / static_psi;
          }
          if (*conformal_state >= 3) {
            psixx[ind] = pxx / static_psi;
            psixy[ind] = pxy / static_psi;
            psixz[ind] = pxz / static_psi;
            psiyy[ind] = pyy / static_psi;
            psiyz[ind] = pyz / static_psi;
            psizz[ind] = pzz / static_psi;
          }

        } /* if conformal-state>0 */
          
        gxx[ind] = pow4 (psi1 / static_psi);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = pow4 (psi1 / static_psi);
        gyz[ind] = 0;
        gzz[ind] = pow4 (psi1 / static_psi);

        kxx[ind] = Aij[0][0] / pow2(psi1);
        kxy[ind] = Aij[0][1] / pow2(psi1);
        kxz[ind] = Aij[0][2] / pow2(psi1);
        kyy[ind] = Aij[1][1] / pow2(psi1);
        kyz[ind] = Aij[1][2] / pow2(psi1);
        kzz[ind] = Aij[2][2] / pow2(psi1);

      }
    }
  }

  free_dvector (F, 0, ntotal - 1);
  free_derivs (&u, ntotal);
  free_derivs (&v, ntotal);
}
