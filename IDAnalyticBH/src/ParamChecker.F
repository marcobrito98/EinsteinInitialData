 /*@@
   @file      ParamChecker.F
   @date      March 1999
   @author    Gabrielle Allen
   @desc 
      Check parameters for black hole initial data and give some
      information
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_arguments.h"
#include "cctk_parameters.h"

#include "CactusEinstein/Einstein/src/Einstein.h"

 /*@@
   @routine    ParamChecker
   @date       March 1999
   @author     Gabrielle Allen
   @desc 
      Check parameters for black hole initial data and give some
      information
   @enddesc 
   @calls     
   @history 
 
   @endhistory 

@@*/

      subroutine ParamChecker(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer CCTK_Equals
      character*80 infoline

      if (CCTK_Equals(initial_data,"schwarzschild") == 1) then

         write(infoline,'(A26,G12.7)') 
     &        'One Black Hole: Throat at ',mass/2d0
         call CCTK_INFO(infoline)

      else if (CCTK_Equals(initial_data,"bl") == 1) then

         write(*,*) "Setting up Brill Lindquist Solution"
         write(*,*) "  Number of black holes (bl_nbh) :",bl_nbh
         write(*,*) "  Black hole masses (bl_M_?) :",bl_M_1,bl_M_2,
     &        bl_M_3,bl_M_4

      else if (CCTK_Equals(initial_data,"multiple_misner") == 1) then

         call CCTK_INFO("Setting up Misner solution for multiple holes")

      else if (CCTK_Equals(initial_data,"misner")==1) then

         call CCTK_INFO("Setting up Misner solution for two holes");
         write(*,*) "  mu = ",mu

      end if

c     Remind users about the conformal metric
c     ---------------------------------------
      if (use_conformal == 1) then
        call CCTK_INFO("Implements conformal metric")
        if (use_conformal_derivs == 1) then
          call CCTK_INFO("  ... and conformal derivatives")
        else
          call CCTK_INFO("  ... but no conformal derivatives")
        end if
      else 
        call CCTK_INFO("Implements non-conformal metric");
        call CCTK_INFO("  (Not usually a good idea!)")        
      end if

      return
      end

