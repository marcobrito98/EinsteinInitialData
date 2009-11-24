#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



extern "C"
void ID_Mag_NS_check_parameters (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (not CCTK_EQUALS (initial_data , "ID_Mag_NS") or
      not CCTK_EQUALS (initial_lapse, "ID_Mag_NS") or
      not CCTK_EQUALS (initial_shift, "ID_Mag_NS") or
      not CCTK_EQUALS (initial_hydro, "ID_Mag_NS") or
      not CCTK_EQUALS (initial_Bmag,  "ID_Mag_NS"))
  {
    CCTK_PARAMWARN ("The parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, HydroBase::initial_hydro, and HydroBase::initial_Bmag must all be set to the value \"ID_Mag_NS\"");
  }
}
