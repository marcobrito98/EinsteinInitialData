
#include "cctk.h"
#include "cctk_parameters.h"
#include "CactusEinstein/Einstein/src/Slicing.h"

void Exact_RegisterSlicing(void) 
{
  int handle;
  handle=Einstein_RegisterSlicing("exact");
}
  
