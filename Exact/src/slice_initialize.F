
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine slice_initialize(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

C     Initialize position of slice in exact solution spacetime.
C     For now the spatial coordinates are identified with those of the 
C     exact solution, and the exact solution time is a Gaussian.

      slicex = x
      slicey = y
      slicez = z
      slicet = slice_gauss_ampl*exp(-(x**2 + y**2 + z**2) 
     $     /slice_gauss_width**2)

C     Calculate Cauchy data and dx^A/dt.

      call slice_data(CCTK_FARGUMENTS)

C     This tells the code that what we have set is the physical metric,
C     not a conformally rescaled one.

      conformal_state = NOCONFORMAL_METRIC
      psi = 1.0D0

      return
      end


