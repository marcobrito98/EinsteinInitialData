#include "cctk.h"
#include "cctk_arguments.h"

#include "Symmetry.h"

void BrilData_InitSymBound(CCTK_CARGUMENTS)
{ 
  DECLARE_CCTK_CARGUMENTS
  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  


  SetCartSymmetry(cctkGH, sym,"IDBrillData::brillpsi");
}
