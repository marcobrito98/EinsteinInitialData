#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine brilldata(CCTK_FARGUMENTS)

c     Author: Carsten Gundlach.
c
c     Driver routine for calculating Brill wave initial data.

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

c     Declare arrays for the linear elliptic solver call:

      integer   Mlin_index,Nsrc_index,field_index,ierr
      integer   metpsi_index(7)

      CCTK_REAL AbsTol(3), RelTol(3)

c     Set up background metric and coefficients for linear solve.

      if (axisym.eq.1) then
         call setupbrilldata2D(CCTK_FARGUMENTS)
      else
         call setupbrilldata3D(CCTK_FARGUMENTS)
      end if

c     Call the Linear Elliptic solver interface
c     to find conformal factor.

      call CCTK_VarIndex (metpsi_index(1), "einstein::gxx")
      call CCTK_VarIndex (metpsi_index(2), "einstein::gxy")
      call CCTK_VarIndex (metpsi_index(3), "einstein::gxz")
      call CCTK_VarIndex (metpsi_index(4), "einstein::gyy")
      call CCTK_VarIndex (metpsi_index(5), "einstein::gyz")
      call CCTK_VarIndex (metpsi_index(6), "einstein::gzz")
      call CCTK_VarIndex (metpsi_index(7), "einstein::psi")

      call CCTK_VarIndex (field_index,     "IDBrillData::brillpsi")
      call CCTK_VarIndex (Mlin_index,      "IDBrillData::Mlinear")
      call CCTK_VarIndex (Nsrc_index,      "IDBrillData::Nsource")

      AbsTol(1)= brill_thresh
      AbsTol(2)= -1
      AbsTol(3)= -1
 
      RelTol(1)= -1
      RelTol(2)= -1
      RelTol(3)= -1

      if (CCTK_EQUALS(brill_bound,"const")) then
         call Ell_SetRealKey(ierr,brill_const_v0,
     .   "EllLinConfMetric::Bnd::Const::V0")
      end if

      if (CCTK_EQUALS(brill_bound,"robin")) then
         call Ell_SetIntKey(ierr,brill_robin_falloff,
     .   "EllLinConfMetric::Bnd::Robin::falloff")
         call Ell_SetRealKey(ierr,brill_robin_inf,
     .   "EllLinConfMetric::Bnd::Robin::inf")
      endif

      if (CCTK_EQUALS(brill_solver,"sor")) then
         call Ell_LinConfMetricSolver(ierr,cctkGH,metpsi_index,
     .   field_index,Mlin_index,Nsrc_index,AbsTol,RelTol,"sor")
      end if

      if (CCTK_EQUALS(brill_solver,"petsc)) then
         call Ell_LinConfMetricSolver(ierr,cctkGH,metpsi_index,
     .   field_index,Mlin_index,Nsrc_index,AbsTol,RelTol,"petsc")
      end if

c     Synchronize conformal factor.

      call CCTK_SyncGroup(cctkGH,"einstein::confac")

c     Reconstruct physical metric.

      call finishbrilldata(CCTK_FARGUMENTS)

      return
      end
