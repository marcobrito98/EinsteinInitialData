/* startup routine for Exact thorn */
/* $Header$ */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "CactusEinstein/CoordGauge/src/Slicing.h"

/*
 * prototypes for scheduled routines
 */
void Exact__RegisterSlicing(CCTK_ARGUMENTS);

/******************************************************************************/

void Exact__RegisterSlicing(CCTK_ARGUMENTS) 
{
  int handle;
  handle=Einstein_RegisterSlicing("exact");
  if (handle<0) CCTK_WARN(1,"Cannot register exact slicing");
}
