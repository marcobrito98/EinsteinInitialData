 /*@@
   @file      ParamChecker.F
   @date      March 1999
   @author    Gabrielle Allen
   @desc 
      Check parameters for black hole initial data and give some
      information
   @enddesc 
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_arguments.h"
#include "cctk_parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

 /*@@
   @routine    ParamChecker
   @date       March 1999
   @author     Gabrielle Allen
   @desc 
      Check parameters for black hole initial data and give some
      information
   @enddesc 
   @calls     
   @history 
 
   @endhistory 

@@*/

void ParamChecker(CCTK_CARGUMENTS)
{
  DECLARE_CCTK_CARGUMENTS
  DECLARE_CCTK_PARAMETERS

  char *message;

  if (CCTK_Equals(initial_data,"schwarzschild") == 1)
  {
    CCTK_INFO("Schwarzschild black hole");
    message = (char *)malloc(200*sizeof(char));
    sprintf(message,"  throat at %f",mass/2.0);
    CCTK_INFO(message);
    free(message);
  }
  else if (CCTK_Equals(initial_data,"bl_bh") == 1)
  {
    CCTK_INFO("Brill Lindquist black holes");
    message = (char *)malloc(200*sizeof(char));
    sprintf(message,"  %d black holes",bl_nbh);
    CCTK_INFO(message);
    if (bl_nbh > 0)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
	      bl_M_1,bl_x0_1,bl_y0_1,bl_z0_1);
      CCTK_INFO(message);
    }
    if (bl_nbh > 1)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
	      bl_M_2,bl_x0_2,bl_y0_2,bl_z0_2);
      CCTK_INFO(message);
    }
    if (bl_nbh > 2)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
	      bl_M_3,bl_x0_3,bl_y0_3,bl_z0_3);
      CCTK_INFO(message);
    }
    if (bl_nbh > 3)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
	      bl_M_4,bl_x0_4,bl_y0_4,bl_z0_4);
      CCTK_INFO(message);
    }
    free(message);
  }
  else if (CCTK_Equals(initial_data,"multiple_misner") == 1)
  {
    CCTK_INFO("Setting up Misner solution for multiple holes");
  }
  else if (CCTK_Equals(initial_data,"misner")==1)
  {
    CCTK_INFO("Two Misner black holes");
    message = (char *)malloc(200*sizeof(char));
    sprintf(message,"  mu is %f",mu);
    CCTK_INFO(message);
    free(message);
  }

  /*     Remind users about the conformal metric
   *     ---------------------------------------
   */
  
  if (use_conformal == 1)
  {
    CCTK_INFO("Black hole initial data uses conformal metric");
    if (use_conformal_derivs == 1)
    {
      CCTK_INFO("  ... and conformal derivatives");
    }
    else
    {
      printf("  ... but no conformal derivatives");
    }
  }
  else
  { 
    printf("Implements non-conformal metric\n");
    printf("  (Not usually a good idea!)\n");
  }        
}

