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

void set_initial_guess(CCTK_POINTER_TO_CONST cctkGH,
                       derivs v)
{
  DECLARE_CCTK_PARAMETERS;

  int nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;
  CCTK_REAL *s_x, *s_y, *s_z;
  CCTK_REAL al, A, Am1, be, B, phi, R, r, X;
  CCTK_INT i, j, k, i3D, ivar, indx;
  derivs U;
  FILE *debug_file;
  CCTK_REAL tmp_r;

  s_x    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_y    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_z    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  allocate_derivs (&U, nvar);
  for (i = 0; i < n1; i++)
    for (j = 0; j < n2; j++)
      for (k = 0; k < n3; k++)
      {
        i3D = Index(0,i,j,k,1,n1,n2,n3);

        al = Pih * (2 * i + 1) / n1;
        A = -cos (al);
        be = Pih * (2 * j + 1) / n2;
        B = -cos (be);
        phi = 2. * Pi * k / n3;

        // Calculation of (X,R)
        AB_To_XR (nvar, A, B, &X, &R, U);
        // Calculation of (x,r)
        C_To_c (nvar, X, R, &(s_x[i3D]), &r, U);
        // Calculation of (y,z)
        rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]), U);
      }
  Set_Initial_Guess_for_u(cctkGH, n1*n2*n3, v.d0, s_x, s_y, s_z);
  for (i = 0; i < n1; i++)
    for (j = 0; j < n2; j++)
      for (k = 0; k < n3; k++)
      {
        indx = Index(0,i,j,k,1,n1,n2,n3);
        v.d0[indx]/=(-cos(Pih * (2 * i + 1) / n1)-1.0);
      }
  Derivatives_AB3 (nvar, n1, n2, n3, v);
  if (do_initial_debug_output)
  {
    debug_file=fopen("initial.dat", "w");
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
      {
        fprintf(debug_file, "%.8g %.8g %.8g %.8g %.8g\n",
                s_x[indx], s_y[indx],
                v.d0[indx],
                (-cos(Pih * (2 * i + 1) / n1)-1.0),
                v.d1[indx]);
      }
    fclose(debug_file);
  }
  free(s_z);
  free(s_y);
  free(s_x);
  free_derivs (&U, nvar);
}

// -------------------------------------------------------------------
void
TwoPunctures (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;

  int antisymmetric_lapse, averaged_lapse;

  int nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int i, j, k, ntotal = n1 * n2 * n3 * nvar;
  static double *F = NULL;
  static derivs u, v;
  double admMass;

  if (! F) {
    /* Solve only when called for the first time */
    F = dvector (0, ntotal - 1);
    allocate_derivs (&u, ntotal);
    allocate_derivs (&v, ntotal);

    CCTK_INFO ("Solving puncture equation");
    /* initialise to 0 */
    for (j = 0; j < ntotal; j++)
      v.d0[j] = 0.0;
    /* call for external initial guess */
    if (use_external_initial_guess)
    {
      set_initial_guess(cctkGH, v);
    }

    Newton (cctkGH, nvar, n1, n2, n3, v, Newton_tol, Newton_maxit);

    F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);
    
    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    admMass = (par_m_plus + par_m_minus
               - 4*par_b*PunctEvalAtArbitPosition(v.d0, 1, 0, 0, n1, n2, n3));
    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass is %g\n", admMass);
  }

  if (CCTK_EQUALS(grid_setup_method, "Taylor expansion"))
  {
    gsm = GSM_Taylor_expansion;
  }
  else if (CCTK_EQUALS(grid_setup_method, "evaluation"))
  {
    gsm = GSM_evaluation;
  }
  else
  {
    CCTK_WARN (0, "internal error");
  }

  antisymmetric_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-antisymmetric");
  averaged_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-averaged");

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
        
        double U;
        switch (gsm)
        {
        case GSM_Taylor_expansion:
          U = PunctTaylorExpandAtArbitPosition
            (0, nvar, n1, n2, n3, v, x[ind], y[ind], z[ind]);
          break;
        case GSM_evaluation:
          U = PunctIntPolAtArbitPosition
            (0, nvar, n1, n2, n3, v, x[ind], y[ind], z[ind]);
          break;
        default:
          assert (0);
        }
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

        } /* if conformal-state > 0 */
          
        puncture_u[ind] = U;

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

        if (antisymmetric_lapse || averaged_lapse) {
/*           const double alp1 = ((1.0 - 0.5 * par_m_plus / r_plus) */
/*                                / (1.0 + 0.5 * par_m_plus / r_plus)); */
/*           const double alp2 = ((1.0 - 0.5 * par_m_minus / r_minus) */
/*                                / (1.0 + 0.5 * par_m_minus / r_minus)); */
/*           alp[ind] = alp1 * alp2; */
          alp[ind] =
            ((1.0 - 0.5 * par_m_plus / r_plus - 0.5 * par_m_minus / r_minus)
             / (1.0 + 0.5 * par_m_plus / r_plus + 0.5 * par_m_minus / r_minus));
          
          if (averaged_lapse) {
            alp[ind] = 0.5 * (1.0 + alp[ind]);
          }
        }

      }
    }
  }

  if (use_sources && rescale_sources)
  {
    Rescale_Sources(cctkGH,
                    cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
                    x, y, z,
                    psi,
                    gxx, gyy, gzz,
                    gxy, gxz, gyz);
  }

  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal - 1);
    free_derivs (&u, ntotal);
    free_derivs (&v, ntotal);
  }
}

