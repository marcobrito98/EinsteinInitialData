 /*@@
   @file      ParamCheck.c
   @date      Ues April 30 2002
   @author    Gabrielle Allen
   @desc 
   Parameter checking stuff for Exact
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(AEIDevelopment_Exact_ParamCheck_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void Exact_ParamCheck(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    Exact_ParamCheck
   @date       Tue April 30 2002
   @author     Tom Goodale
   @desc 
   Scheduled routine to detect invalid parameter settings.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @hdate      Tue Jun 11 18:25:49 CEST 2002
   @hauthor    Jonathan Thornburg
   @desc       Add a test that we're not trying to set the shift
               with storage *not* present for it
   @enddesc
   @endhistory 
 
 @@*/
void Exact_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  if(! CCTK_EQUALS(metric_type, "physical"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\"");
  }

  if ( CCTK_EQUALS(shift_evolution_method, "exact")
       && CCTK_EQUALS(initial_shift, "none") )
  {
    CCTK_PARAMWARN("can't set the shift if there's no storage for it!");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
