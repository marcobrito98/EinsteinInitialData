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
  if( ! CCTK_Equals(metric_type, "physical") && 
      ! CCTK_Equals(metric_type, "static conformal"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\" and \"static conformal\"");
  }

  if (CCTK_Equals(metric_type, "static conformal"))
  {
    if (!CCTK_Equals(conformal_storage,"factor"))
    {
      CCTK_PARAMWARN("BrillData only sets the conformal factor (this could easily be changed ... please ask)");
    }
  }

  CCTK_INFO("Setting up Brill data");
  if (CCTK_Equals(metric_type, "static conformal"))
  {
     CCTK_VInfo(CCTK_THORNSTRING,"  ... using trivial conformal %s",conformal_storage);
  }
  else if (CCTK_Equals(metric_type, "physical"))
  {
    CCTK_INFO("  ... using physical metric");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
