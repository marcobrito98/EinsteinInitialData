
#include "cctk.h"
#include "cctk_arguments.h"
#include "CactusEinstein/Einstein/src/Einstein.h"

      subroutine slice_initialize(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS

      CCTK_REAL slice_gauss_ampl, slice_gauss_width

C     Initialize position of slice in exact solution spacetime.
C     For now something simple.

      slicex = x
      slicey = y
      slicez = z
      slicet = slice_gauss_ampl*exp(-(x**2 + y**2 + z**2) 
     $     /slice_gauss_width**2)

C     Calculate Cauchy data and dx^A/dt.

      call slice_data(CCTK_FARGUMENTS)

C     Extrapolate g, K, etc. to the boundary, where slice_data cannot
C     calculate them, because it cannot take derivatives there. At all
C     other timesteps cactus will do this automatically.

      call boundaries(CCTK_FARGUMENTS)

C     This wasted a day of work. PW found it...

      conformal_state = NOCONFORMAL_METRIC
      psi = 1.0D0

      return
      end


