#include "cctk.h"
#include "declare_arguments.h"
#include "declare_parameters.h"

      subroutine ParamChecker(CCTK_FARGUMENTS)

      implicit none

      DECLARE_CCTK_FARGUMENTS
      DECLARE_PARAMETERS

      integer CCTK_Equals

      if (CCTK_Equals(initial_data,"schwarzschild")==1) then
         write(*,*) 'One Black Hole:  Throat at :',mass/2d0
         if (use_conformal == 1) then
            write(*,*) "Uses conformal metric"
         else 
            write(*,*) "Uses non-conformal metric"         
         end if
      else if (CCTK_Equals(initial_data,"bl")==1) then
         write(*,*) "Setting up Brill Lindquist Solution"
         write(*,*) "   Uses conformal metric"
         write(*,*) "   Number of black holes (bl_nbh) :",bl_nbh
         write(*,*) "   Black hole masses (bl_M_?) :",bl_M_1,bl_M_2,
     &        bl_M_3,bl_M_4
      else if (CCTK_Equals(initial_data,"misner")==1) then
         write(*,*) "Setting up Misner Solution"
         write(*,*) "   Uses conformal metric"
      end if

      return
      end

