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
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

#define SQR(a) ((a)*(a))

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
   @calls     
   @history 
 
   @endhistory 

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
void Misner_standard(CCTK_CARGUMENTS)
{
  DECLARE_CCTK_CARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i,j,k;
  int n;
  char *message;
  CCTK_REAL zero, one, three;
  CCTK_REAL csch, coth, r1, r2;
  CCTK_REAL x_squared, y_squared;
  CCTK_REAL r1_cubed, r2_cubed;
  CCTK_REAL r1_5, r2_5;
  CCTK_INT powfac;

  int nx,ny,nz;
  int index;
  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];


  zero  = 0.0;
  one   = 1.0;
  three = 3.0;

  /*     Initialize so we can accumulate
   *     ------------------------------- 
   */
  for(k=0; k < nz; k++)
  {
    for(j=0; j < ny; j++)
    {
      for(i=0; i < nx; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH, i,j,k);
        
        psi[index] = one;
      }
    }
  }

  if (use_conformal_derivs == 1)
  {
    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          
          psix[index]  = zero;
          psiy[index]  = zero;
          psiz[index]  = zero;
          psixx[index] = zero;
          psixy[index] = zero;
          psixz[index] = zero;
          psiyy[index] = zero;
          psiyz[index] = zero;
          psizz[index] = zero;
        }
      }
    }
  }

  for(n = nmax; n >= 1; n--)
  {
    csch = one/sinh(mu*n);
    coth = one/tanh(mu*n);
    
    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          
          x_squared = SQR(x[index]);
          y_squared = SQR(y[index]);
          r1 = sqrt(x_squared+y_squared+SQR(z[index]+coth));
          r2 = sqrt(x_squared+y_squared+SQR(z[index]-coth));
          
          psi[index] += csch*(one/r1 + one/r2);
          
          if (use_conformal_derivs == 1)
          {
            r1_cubed = r1*r1*r1;
            r2_cubed = r2*r2*r2;
            r1_5     = pow(r1, 5);
            r2_5     = pow(r2, 5);
            psix[index]  +=  (-(x[index]/r2_cubed)-x[index]/r1_cubed)*csch;
            psiy[index]  +=  (-(y[index]/r2_cubed)-y[index]/r1_cubed)*csch;
            psiz[index]  +=  (-((z[index]-coth)/r2_cubed)-(z[index] + coth)/r1_cubed)*csch;
            psixx[index] +=  ((three*x_squared)/r2_5-one/r2_cubed
                              +(three*x_squared)/r1_5-one/r1_cubed)*csch;
            psixy[index] +=  ((three*x[index]*y[index])/r2_5
                              +(three*x[index]*y[index])/r1_5)*csch;
            psixz[index] +=  ((three*x[index]*(z[index]-coth))/r2_5 
                              +(three*x[index]*(z[index]+coth))/r1_5)*csch;
            psiyy[index] +=  ((three*y_squared)/r2_5-one/r2_cubed 
                              +(three*y_squared)/r1_5-one/r1_cubed)*csch;
            psiyz[index] +=  ((three*y[index]*(z[index]-coth))/r2_5
                              +(three*y[index]*(z[index]+coth))/r1_5)*csch;
            psizz[index] += (-one/r2_cubed+(three*SQR(z[index]-coth))/r2_5+
                             (three*SQR(z[index]+coth))/r1_5-one/r1_cubed)*csch;
          }
        }
      }
    }
  }

  /*     Cactus convention
   *     -----------------
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
  
  /*     compute the ADM mass
   *     -------------------- 
   */
  mass = zero;
  
  for(n = 1; n <= nmax; n++)
  {
    mass += 4./sinh(n*mu);
  }

  message = (char *)malloc(200*sizeof(char));
  sprintf(message,"ADM mass is %f",mass);
  CCTK_INFO(message);
  free(message);
  
  /*     Should initialize lapse to Cadez value if possible
   *     --------------------------------------------------
   */

  if (CCTK_Equals(initial_lapse,"cadez") )
  { 
    CCTK_INFO("Initialise with cadez lapse");
    
    for(k=0; k < nz; k++)
    {
      for(j=0; j < ny; j++)
      {
        for(i=0; i < nx; i++)
        {
          index = CCTK_GFINDEX3D(cctkGH, i,j,k);
          
          x_squared = SQR(x[index]);
          y_squared = SQR(y[index]);
            
          alp[index] = one;

          powfac = 1;

          for(n = 1; n <= nmax; n++)
          {
            coth = one/tanh(mu*n);
            r1 = sqrt(x_squared+y_squared+SQR(z[index]+coth));
            r2 = sqrt(x_squared+y_squared+SQR(z[index]-coth));
            powfac *= -1;

            alp[index] += powfac * one/sinh(mu*n)*(one/r1 + one/r2);
          }
          
          alp[index] /= psi[index];
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

  USE_CCTK_PARAMETERS
  USE_CCTK_CARGUMENTS

}
