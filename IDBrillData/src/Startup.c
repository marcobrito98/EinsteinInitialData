#include "cctk.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

void BrilData_InitSymBound(CCTK_ARGUMENTS);

void BrilData_InitSymBound(CCTK_ARGUMENTS)
{ 
  DECLARE_CCTK_ARGUMENTS
  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  
  SetCartSymVN(cctkGH, sym,"IDBrillData::brillpsi");

  return;
} 
