// TwoPunctures:  File  "TwoPunctures.c"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

// -------------------------------------------------------------------
void
TwoPunctures_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (par_use_sources)
  {
    CCTK_INFO("Solving for BH-NS");
    if (CCTK_IsFunctionAliased ("Set_Rho_ADM"))
      CCTK_INFO("Aliased Functions found");
    else
      CCTK_WARN(0, "I found no (aliased) function for matter sources, but "
                   "was said to use matter.\n");
  }
  else
    CCTK_INFO("not using sources (only BHs)");
}
