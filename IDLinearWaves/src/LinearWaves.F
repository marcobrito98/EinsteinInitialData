c/*@@
c  @file      linearWaves.F
c  @date      
c  @author    Joan Masso
c  @desc 
c             Driver for the analytical (or almost) linear wave models 
c  @enddesc 
c@@*/

c/*@@
c  @routine    linearWaves
c  @date       
c  @author     Joan Masso 
c  @desc 
c  
c  @enddesc 
c  @calls     planewaves teukwaves
c  @calledby   
c  @history Converted to Cactus4.0. the code and its comments of the
   original authors are kept. There is NO BM support at this point.
   Id like to know how Cactus handles this w/o #ifdef !
c
c  @endhistory 
c@@*/

#include "cctk.h"
#include "declare_arguments.h"
#include "declare_parameters.h"

      subroutine linearWaves(CCTK_FARGUMENTS)

      implicit none
      
      DECLARE_CCTK_FARGUMENTS
      DECLARE_PARAMETERS

      integer CCTK_Equals

c     Call the Teukolsky wave init routine
      if (CCTK_Equals(initial_data,'teukwaves').eq.1) then

         call teukwaves(CCTK_FARGUMENTS)

c     at this point CACTUS 3.2.0 has the following calls
c     #ifdef BMUTIL call metricderiv
c     #ifdef BONAMASSO: vxyz=0; call vectorini;

      elseif (CCTK_Equals(initial_data,'planewaves') == 1) then

         call planewaves(CCTK_FARGUMENTS)

c     at this point CACTUS 3.2.0 has the following calls     
c     #ifdef BMUTIL call vectorini;

      else 

         call CCTK_Warn(0,"Thorn linearWaves: model initialization failed");
      endif
      
      return
      end
