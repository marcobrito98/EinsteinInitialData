
#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "CactusEinstein/Einstein/src/Einstein.h"
#include "CactusElliptic/EllBase/src/EllBase.h"

      subroutine brilldata(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer ipsi,iMcoeff,iNcoeff
      integer ierr

      CCTK_REAL AbsTol(3),RelTol(3)

c     Get indices for grid functions.

      call CCTK_VarIndex(ipsi   , "idbrilldata::brillpsi")
      if (ipsi .lt. 0) then
	call CCTK_WARN(0,"Grid variable index for iphi not found")
      end if

      call CCTK_VarIndex(iMcoeff, "idbrilldata::brillMlinear")
      if (iMcoeff .lt. 0) then
	call CCTK_WARN(0,"Grid variable index for Mcoeff not found")
      end if

      call CCTK_VarIndex(iNcoeff, "idbrilldata::brillNsource")
      if (iNcoeff .lt. 0) then
	call CCTK_WARN(0,"Grid variable index for Ncoeff not found")
      end if

c     Set up background metric and coefficients for linear solve.

      if (axisym.eq.1) then
         call setupbrilldata2D(CCTK_FARGUMENTS)
      else
         call setupbrilldata3D(CCTK_FARGUMENTS)
      end if

c     Tolerances for elliptic solve.

      AbsTol(1)= brill_thresh
      AbsTol(2)= brill_thresh
      AbsTol(3)= brill_thresh
 
      RelTol(1)= -1
      RelTol(2)= -1
      RelTol(3)= -1

c     Elliptic solver.

      if (axisym.eq.1) then

         call Ell_SetRealKey(ierr, 1.0d0, "EllLinFlat::Bnd::Robin::inf")
         call Ell_SetIntKey (ierr, 1,     "EllLinFlat::Bnd::Robin::falloff")

         if (CCTK_EQUALS(brill_solver,"sor")) then
            call Ell_LinFlatSolver(ierr,cctkGH,ipsi,
     .      iMcoeff,iNcoeff,AbsTol,RelTol,"sor")
         end if

         if (CCTK_EQUALS(brill_solver,"petsc")) then
            call Ell_LinFlatSolver(ierr,cctkGH,ipsi,
     .      iMcoeff,iNcoeff,AbsTol,RelTol,"petsc")
         end if

         if (CCTK_EQUALS(brill_solver,"bam")) then
            call Ell_LinFlatSolver(ierr,cctkGH,ipsi,
     .      iMcoeff,iNcoeff,AbsTol,RelTol,"bam")
         end if

      else

      end if

c     Check for errors.

      if (ierr .eq. ELL_SUCCESS) then
        call CCTK_INFO("Leaving elliptic solver: solve successful")
      else if (ierr .eq. ELL_NOCONVERGENCE) then
	call CCTK_INFO("Leaving elliptic solver: solver failed to converge")
      else if (ierr .eq. ELL_NOSOLVER) then
	call CCTK_INFO("Elliptic solver not found")
      else 
        call CCTK_INFO("Leaving elliptic solver: solve failed")	
      end if

c     Synchronization and symmetry boundaries.

      call CCTK_SyncGroup(cctkGH,"idbrilldata::brillconf")
      call CartSymGN(ierr,cctkGH,"idbrilldata::brillconf")

c     Reconstruct physical metric.

      call finishbrilldata(CCTK_FARGUMENTS)

      return
      end
