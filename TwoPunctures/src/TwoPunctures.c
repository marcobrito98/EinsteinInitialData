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

  CCTK_INFO ("Beginning elliptic solving");
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
        const double psi1 = 1 + 0.5 * par_m_plus / r_plus
                              + 0.5 * par_m_minus / r_minus + U;
        const double psi4 = pow2(pow2(psi1));
        const double psim2 = 1.0/pow2(psi1);

        double Aij[3][3];
        BY_Aijofxyz (x[ind], y[ind], z[ind], Aij);

        switch (*conformal_state) {
        case 0:
          gxx[ind] = psi4;
          gxy[ind] = 0;
          gxz[ind] = 0;
          gyy[ind] = psi4;
          gyz[ind] = 0;
          gzz[ind] = psi4;
          break;
          
        case 3:
          /* not yet supported */
          assert (0);
          /* fall through */

        case 2:
          /* not yet supported */
          assert (0);
          /* fall through */

        case 1:
          psi[ind] = psi1;

          gxx[ind] = 1;
          gxy[ind] = 0;
          gxz[ind] = 0;
          gyy[ind] = 1;
          gyz[ind] = 0;
          gzz[ind] = 1;
          break;

        default:
          assert(0);
        }

        kxx[ind] = psim2 * Aij[0][0];
        kxy[ind] = psim2 * Aij[0][1];
        kxz[ind] = psim2 * Aij[0][2];
        kyy[ind] = psim2 * Aij[1][1];
        kyz[ind] = psim2 * Aij[1][2];
        kzz[ind] = psim2 * Aij[2][2];

      }
    }
  }

  free_dvector (F, 0, ntotal - 1);
  free_derivs (&u, ntotal);
  free_derivs (&v, ntotal);
}
