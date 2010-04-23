 /*@@
   @file      Hydro_InitExcisionMask.c
   @date      Wed Mar 13 14:18:38 2002
   @author    Andrea Nerozzi, Frank Loeffler
   @desc 
   Creates the mask for an excision region 
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include <stdio.h>
#include <math.h>

 /*@@
   @routine    Excision Mask 
   @date       Mon Apr 7 14:18:38 2003
   @author     Andrea Nerozzi
   @desc 
   Sets up the excision region
   @enddesc 
   @calls     
   @calledby   
   @history    added 'sphere' and rewrote to be multi-processor-safe
 
   @endhistory 

@@*/

void Hydro_InitExcisionMask(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT  i,j,k, nx, ny, nz, point;
  CCTK_INT excised, normal;
  CCTK_REAL emask_excised, emask_normal;
  CCTK_REAL x_min, x_max, y_min, y_max, z_min, z_max;
  CCTK_REAL x_size, y_size, z_size;
  CCTK_REAL x_frac, y_frac, z_frac;

  if (hydro_initexcision)
  {
    CCTK_INFO("Setting up Hydro Initial Excision region.");
  }

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];

  CCTK_CoordRange(cctkGH, &x_min, &x_max, 1, NULL, "cart3d");
  CCTK_CoordRange(cctkGH, &y_min, &y_max, 2, NULL, "cart3d");
  CCTK_CoordRange(cctkGH, &z_min, &z_max, 3, NULL, "cart3d");
  x_size = x_max - x_min;
  y_size = y_max - y_min;
  z_size = z_max - z_min;

  excised=1;
  normal =0;
  emask_excised=0.0;
  emask_normal=1.0;
  if (hydro_initexcision_flip)
  {
    excised=0;
    normal =1;
    emask_excised=1.0;
    emask_normal=0.0;
  }

  for (k=0 ; k<nz ; k++)
  {
    for (j=0 ; j<ny ; j++)
    {
      for (i=0 ; i<nx; i++)
      { 
        point = CCTK_GFINDEX3D(cctkGH,i,j,k);
        x_frac = (x[point]-x_min) / x_size + 1.e-10;
        y_frac = (y[point]-y_min) / y_size + 1.e-10;
        z_frac = (z[point]-z_min) / z_size + 1.e-10;

        if (hydro_initexcision)
        {
          if (hydro_initexcision_old_mask)
            emask[point]=emask_normal;
          if (CCTK_Equals(hydro_initexcision_type,"box"))
          {
              if (CCTK_Equals(domain,"full"))
              {
                 if ( (hydro_initexcision_coordinate_length <= 0.0) &&
                      ( ( x_frac > 0.5 - hydro_initexcision_fraction) &&
                        ( x_frac < 0.5 + hydro_initexcision_fraction) &&
                        ( y_frac > 0.5 - hydro_initexcision_fraction) &&
                        ( y_frac < 0.5 + hydro_initexcision_fraction) &&
                        ( z_frac > 0.5 - hydro_initexcision_fraction) &&
                        ( z_frac < 0.5 + hydro_initexcision_fraction)
                      ) ||
                      ( (hydro_initexcision_coordinate_length > 0.0) &&
                        ( fabs(x[point]-hydro_initexcision_position_x) <=
                          hydro_initexcision_coordinate_length*0.5) &&
                        ( fabs(y[point]-hydro_initexcision_position_y) <=
                          hydro_initexcision_coordinate_length*0.5) &&
                        ( fabs(z[point]-hydro_initexcision_position_z) <=
                          hydro_initexcision_coordinate_length*0.5)
                      )
                    )
                 {
                    hydro_excision_mask[point] = excised;
                    if (hydro_initexcision_old_mask)
                      emask[point]=emask_excised;
                 }
                 else
                 {
                    hydro_excision_mask[point] = normal;
                 }
              }
              else if (CCTK_Equals(domain,"octant"))
              {
                 if ( ( (hydro_initexcision_coordinate_length <= 0.0) &&
                        ( x_frac < hydro_initexcision_fraction) &&
                        ( y_frac < hydro_initexcision_fraction) &&
                        ( z_frac < hydro_initexcision_fraction)
                      ) || 
                      ( (hydro_initexcision_coordinate_length > 0.0) &&
                        ( fabs(x[point]-hydro_initexcision_position_x) <=
                          hydro_initexcision_coordinate_length*0.5) &&
                        ( fabs(y[point]-hydro_initexcision_position_y) <=
                          hydro_initexcision_coordinate_length*0.5) &&
                        ( fabs(z[point]-hydro_initexcision_position_z) <=
                          hydro_initexcision_coordinate_length*0.5)
                      )
                    )
                 {
                    hydro_excision_mask[point] = excised;
                    if (hydro_initexcision_old_mask)
                      emask[point]=emask_excised;
                 }
                 else
                 {
                    hydro_excision_mask[point] = normal;
                 }
              }
          }
          else if (CCTK_Equals(hydro_initexcision_type,"x-axis"))
          {
            if ( x_frac < 1.0-hydro_initexcision_fraction )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (hydro_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(hydro_initexcision_type,"y-axis"))
          {
            if ( y_frac < 1.0-hydro_initexcision_fraction )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (hydro_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(hydro_initexcision_type,"z-axis"))
          {
            if ( z_frac < 1.0-hydro_initexcision_fraction )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (hydro_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(hydro_initexcision_type,"diagonal"))
          {
            if ( x_frac + y_frac + z_frac <
                 3.0*(1.0-hydro_initexcision_fraction) )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (hydro_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(hydro_initexcision_type,"sphere"))
          {
            if ( (hydro_initexcision_coordinate_length <= 0.0) &&
                 ( ( (x[point]-hydro_initexcision_position_x) *
                     (x[point]-hydro_initexcision_position_x) +
                     (y[point]-hydro_initexcision_position_y) *
                     (y[point]-hydro_initexcision_position_y) +
                     (z[point]-hydro_initexcision_position_z) *
                     (z[point]-hydro_initexcision_position_z) ) * 4.0
                   > (hydro_initexcision_fraction * x_size *
                      hydro_initexcision_fraction * x_size) ) ||
                 ( (hydro_initexcision_coordinate_length > 0) &&
                   ( (hydro_initexcision_coordinate_length < 
                      cctkGH->cctk_delta_space[0]*
                      hydro_initexcision_min_points) ||
                       ( (x[point]-hydro_initexcision_position_x) *
                         (x[point]-hydro_initexcision_position_x) +
                         (y[point]-hydro_initexcision_position_y) *
                         (y[point]-hydro_initexcision_position_y) +
                         (z[point]-hydro_initexcision_position_z) *
                         (z[point]-hydro_initexcision_position_z)
                     > hydro_initexcision_coordinate_length*
                       hydro_initexcision_coordinate_length) ) ) )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (hydro_initexcision_old_mask)
              {
                emask[point]=emask_excised;
              }
            }
          }
        }
        else
        {
          hydro_excision_mask[point] = normal;
          if (hydro_initexcision_old_mask)
            emask[point]=emask_normal;
        }
      }
    }
  }
  return;
}

void Hydro_InitExcisionMask_poststep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (hydro_initexcision_poststep)
      Hydro_InitExcisionMask(CCTK_PASS_CTOC);
}

