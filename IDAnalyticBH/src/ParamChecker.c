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

  if (CCTK_Equals(initial_data,"schwarzschild") == 1)
  {
    printf("One Black Hole:  Throat at : %f\n",mass/2.0);
  }
  else if (CCTK_Equals(initial_data,"bl") == 1)
  {
    printf("Setting up Brill Lindquist Solution\n");
    printf("  Number of black holes (bl_nbh) : %d\n",bl_nbh);
    printf("  Black hole masses (bl_M_?) : %f %f %f %f\n",bl_M_1,bl_M_2,bl_M_3,bl_M_4);
  }
  else if (CCTK_Equals(initial_data,"multiple_misner") == 1)
  {
    printf("Setting up Misner solution for multiple holes\n");
  }
  else if (CCTK_Equals(initial_data,"misner")==1)
  {
    printf("Setting up Misner solution for two holes\n");
    printf("  mu = %f\n");
  }

  /*     Remind users about the conformal metric
   *     ---------------------------------------
   */
  
  if (use_conformal == 1)
  {
    printf("Implements conformal metric\n");
    if (use_conformal_derivs == 1)
    {
      printf("  ... and conformal derivatives\n");
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

