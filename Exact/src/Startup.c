
#include "cctk.h"
#include "cctk_Parameters.h"
#include "CactusEinstein2/CoordGauge/src/Slicing.h"

void Exact_RegisterSlicing(void) 
{
  int handle;
  handle=Einstein_RegisterSlicing("exact");
  if (handle<0) CCTK_WARN(1,"Cannot register exact slicing");
}
