/* $Header$ */

#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static void callback (int idx, const char * optstring, void * callback_arg);

/** Ensure that all ADMBase initial data that are supposed to be read
    from a file are actually scheduled for the file reader.  */
void IDFileADM_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  char * variable_is_read;
  int i;
  int nvars;
  
  variable_is_read = malloc (CCTK_NumVars());
  assert (variable_is_read);
  for (i=0; i<CCTK_NumVars(); ++i) {
    variable_is_read[i] = 0;
  }
  
  nvars = CCTK_TraverseString
    (filereader_ID_files, callback, variable_is_read, CCTK_GROUP_OR_VAR);
  assert (nvars >= 0);
  
  if (CCTK_EQUALS(initial_lapse, "read from file")) {
    int const ialp = CCTK_VarIndex ("ADMBase::alp");
    assert (ialp >= 0);
    if (! variable_is_read[ialp]) {
      CCTK_PARAMWARN ("The lapse is initialised using the file reader, but the group ADMBase::lapse has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_files\" accordingly.");
    }
  }
  
  if (CCTK_EQUALS(initial_shift, "read from file")) {
    int const ibetax = CCTK_VarIndex ("ADMBase::betax");
    int const ibetay = CCTK_VarIndex ("ADMBase::betay");
    int const ibetaz = CCTK_VarIndex ("ADMBase::betaz");
    assert (ibetax >= 0);
    assert (ibetay >= 0);
    assert (ibetaz >= 0);
    if (! variable_is_read[ibetax]
        || ! variable_is_read[ibetay]
        || ! variable_is_read[ibetaz]) {
      CCTK_PARAMWARN ("The shift is initialised using the file reader, but the group ADMBase::shift has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_files\" accordingly.");
    }
  }
  
  if (CCTK_EQUALS(initial_shift, "read from file")) {
    int const igxx = CCTK_VarIndex ("ADMBase::gxx");
    int const igxy = CCTK_VarIndex ("ADMBase::gxy");
    int const igxz = CCTK_VarIndex ("ADMBase::gxz");
    int const igyy = CCTK_VarIndex ("ADMBase::gyy");
    int const igyz = CCTK_VarIndex ("ADMBase::gyz");
    int const igzz = CCTK_VarIndex ("ADMBase::gzz");
    assert (igxx >= 0);
    assert (igxy >= 0);
    assert (igxz >= 0);
    assert (igyy >= 0);
    assert (igyz >= 0);
    assert (igzz >= 0);
    if (! variable_is_read[igxx]
        || ! variable_is_read[igxy]
        || ! variable_is_read[igxz]
        || ! variable_is_read[igyy]
        || ! variable_is_read[igyz]
        || ! variable_is_read[igzz]) {
      CCTK_PARAMWARN ("The metric is initialised using the file reader, but the group ADMBase::metric has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_files\" accordingly.");
    }
  }
  
  free (variable_is_read);
}

/** Mark a variable as to be read from the file reader.  */
static void callback (int idx, const char * optstring, void * callback_arg)
{
  assert (idx>=0 && idx<CCTK_NumVars());
  ((char *)callback_arg)[idx] = 1;
}
