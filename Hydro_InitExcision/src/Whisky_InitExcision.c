 /*@@
   @file      Whisky_ExcisionMask.c
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

void Whisky_InitExcisionMask(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT  i,j,k, nx, ny, nz, point;
  CCTK_INT excised, normal;
  CCTK_REAL emask_excised, emask_normal;
  CCTK_REAL x_min, x_max, y_min, y_max, z_min, z_max;
  CCTK_REAL x_size, y_size, z_size;
  CCTK_REAL x_frac, y_frac, z_frac;

  if (whisky_initexcision)
  {
    CCTK_INFO("Setting up Whisky Initial Excision region.");
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
  if (whisky_initexcision_flip)
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

        if (whisky_initexcision)
        {
          if (whisky_initexcision_old_mask)
            emask[point]=emask_normal;
          if (CCTK_Equals(whisky_initexcision_type,"box"))
          {
              if (CCTK_Equals(domain,"full"))
              {
                 if ( (whisky_initexcision_coordinate_length <= 0.0) &&
                      ( ( x_frac > 0.5 - whisky_initexcision_fraction) &&
                        ( x_frac < 0.5 + whisky_initexcision_fraction) &&
                        ( y_frac > 0.5 - whisky_initexcision_fraction) &&
                        ( y_frac < 0.5 + whisky_initexcision_fraction) &&
                        ( z_frac > 0.5 - whisky_initexcision_fraction) &&
                        ( z_frac < 0.5 + whisky_initexcision_fraction)
                      ) ||
                      ( (whisky_initexcision_coordinate_length > 0.0) &&
                        ( fabs(x[point]-whisky_initexcision_position_x) <=
                          whisky_initexcision_coordinate_length*0.5) &&
                        ( fabs(y[point]-whisky_initexcision_position_y) <=
                          whisky_initexcision_coordinate_length*0.5) &&
                        ( fabs(z[point]-whisky_initexcision_position_z) <=
                          whisky_initexcision_coordinate_length*0.5)
                      )
                    )
                 {
                    hydro_excision_mask[point] = excised;
                    if (whisky_initexcision_old_mask)
                      emask[point]=emask_excised;
                 }
                 else
                 {
                    hydro_excision_mask[point] = normal;
                 }
              }
              else if (CCTK_Equals(domain,"octant"))
              {
                 if ( ( (whisky_initexcision_coordinate_length <= 0.0) &&
                        ( x_frac < whisky_initexcision_fraction) &&
                        ( y_frac < whisky_initexcision_fraction) &&
                        ( z_frac < whisky_initexcision_fraction)
                      ) || 
                      ( (whisky_initexcision_coordinate_length > 0.0) &&
                        ( fabs(x[point]-whisky_initexcision_position_x) <=
                          whisky_initexcision_coordinate_length*0.5) &&
                        ( fabs(y[point]-whisky_initexcision_position_y) <=
                          whisky_initexcision_coordinate_length*0.5) &&
                        ( fabs(z[point]-whisky_initexcision_position_z) <=
                          whisky_initexcision_coordinate_length*0.5)
                      )
                    )
                 {
                    hydro_excision_mask[point] = excised;
                    if (whisky_initexcision_old_mask)
                      emask[point]=emask_excised;
                 }
                 else
                 {
                    hydro_excision_mask[point] = normal;
                 }
              }
          }
          else if (CCTK_Equals(whisky_initexcision_type,"x-axis"))
          {
            if ( x_frac < 1.0-whisky_initexcision_fraction )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (whisky_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(whisky_initexcision_type,"y-axis"))
          {
            if ( y_frac < 1.0-whisky_initexcision_fraction )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (whisky_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(whisky_initexcision_type,"z-axis"))
          {
            if ( z_frac < 1.0-whisky_initexcision_fraction )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (whisky_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(whisky_initexcision_type,"diagonal"))
          {
            if ( x_frac + y_frac + z_frac <
                 3.0*(1.0-whisky_initexcision_fraction) )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (whisky_initexcision_old_mask)
                emask[point]=emask_excised;
            }
          }
          else if (CCTK_Equals(whisky_initexcision_type,"sphere"))
          {
            if ( (whisky_initexcision_coordinate_length <= 0.0) &&
                 ( ( (x[point]-whisky_initexcision_position_x) *
                     (x[point]-whisky_initexcision_position_x) +
                     (y[point]-whisky_initexcision_position_y) *
                     (y[point]-whisky_initexcision_position_y) +
                     (z[point]-whisky_initexcision_position_z) *
                     (z[point]-whisky_initexcision_position_z) ) * 4.0
                   > (whisky_initexcision_fraction * x_size *
                      whisky_initexcision_fraction * x_size) ) ||
                 ( (whisky_initexcision_coordinate_length > 0) &&
                   ( (whisky_initexcision_coordinate_length < 
                      cctkGH->cctk_delta_space[0]*
                      whisky_initexcision_min_points) ||
                       ( (x[point]-whisky_initexcision_position_x) *
                         (x[point]-whisky_initexcision_position_x) +
                         (y[point]-whisky_initexcision_position_y) *
                         (y[point]-whisky_initexcision_position_y) +
                         (z[point]-whisky_initexcision_position_z) *
                         (z[point]-whisky_initexcision_position_z)
                     > whisky_initexcision_coordinate_length*
                       whisky_initexcision_coordinate_length) ) ) )
            {
              hydro_excision_mask[point] = normal;
            }
            else 
            {
              hydro_excision_mask[point] = excised;
              if (whisky_initexcision_old_mask)
              {
                emask[point]=emask_excised;
              }
            }
          }
        }
        else
        {
          hydro_excision_mask[point] = normal;
          if (whisky_initexcision_old_mask)
            emask[point]=emask_normal;
        }
      }
    }
  }
  return;
}

void Whisky_InitExcisionMask_poststep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (whisky_initexcision_poststep)
      Whisky_InitExcisionMask(CCTK_PASS_CTOC);
}

