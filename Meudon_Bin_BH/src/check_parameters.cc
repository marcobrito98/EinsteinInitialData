#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



extern "C"
void ID_Bin_BH_check_parameters (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (not CCTK_EQUALS (initial_data , "ID_Bin_BH") or
      not CCTK_EQUALS (initial_lapse, "ID_Bin_BH") or
      not CCTK_EQUALS (initial_shift, "ID_Bin_BH"))
  {
    CCTK_PARAMWARN ("The parameters ADMBase::initial_data, ADMBase::initial_lapse, and ADMBase::initial_shift must all be set to the value \"ID_Bin_BH\"");
  }
}
