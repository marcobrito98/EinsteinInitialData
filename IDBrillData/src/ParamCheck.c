 /*@@
   @file      ParamCheck.c
   @date      Sun May  5 19:38:51 CEST 2002
   @author    David Rideout
   @desc 
   Check the parameters for IDBrillData
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_IDBrillData_ParamCheck_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void IDBrillData_ParamChecker(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

void IDBrillData_ParamChecker(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  /* Do we know how to deal with this type of metric ? */
  if( ! CCTK_EQUALS(metric_type, "static conformal"))
    if (CCTK_EQUALS(metric_type, "physical"))
      {
	/* This thorn seems to 'return' the physical metric regardless, but will set the conformal factor to 1 if metric_type is static conformal 
	 Is this the desired behavior? */
      } else 
	{
	  CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\" and \"static conformal\"");
	}
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
