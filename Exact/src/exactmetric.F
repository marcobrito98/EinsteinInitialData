#include "cctk.h"
#include "cctk_parameters.h"
#include "cctk_arguments.h"

C     This is just a wrapper for different routines in place of exactmetric.

      subroutine exactmetric(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guxz)

      implicit none

      DECLARE_CCTK_PARAMETERS

      CCTK_REAL x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guxz

      logical firstcall

      integer CCTK_Equals
      
C     Call the corresponding routine.

      if (CCTK_Equals(exactmodel,"boostrot").eq.1) then

         call boostrotmetric(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"finkelstein").eq.1) then

         call finkelstein(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"kerrschild").eq.1) then

         call kerrschild(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"flatschwarz").eq.1) then

         call flatschwarz(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"starschwarz").eq.1) then

         call starschwarz(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"novikov").eq.1) then

         call novikov(x,y,z,t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"bowl").eq.1) then

         call bowl(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"minkowski").eq.1) then

         call minkowski(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"fakebinary").eq.1) then

         call fakebinary(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      elseif (CCTK_Equals(exactmodel,"multiBH").eq.1) then

         call multi_BHs(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz, 
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz, 
     $        guxx, guyy, guzz, guxy, guyz, guxz)

      else

         call CCTK_WARN(0,"Unknown value of parameter exactmodel")

      end if

      return
      end
