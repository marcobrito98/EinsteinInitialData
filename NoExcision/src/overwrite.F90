! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NoExcision_Overwrite (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL :: cx, cy, cz, rad
  integer   :: n
  
  do n = 1, num_regions
     
     cx = centre_x(n)
     cy = centre_y(n)
     cz = centre_z(n)
     rad = radius(n)
     
     if (overwrite_geometry /= 0) then
        
        if (conformal_state >= 1) then
           where ((x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= rad**2)
              psi = 1
           end where
        end if
        if (conformal_state >= 2) then
           where ((x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= rad**2)
              psix = 0
              psiy = 0
              psiz = 0
           end where
        end if
        if (conformal_state >= 3) then
           where ((x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= rad**2)
              psixx = 0
              psixy = 0
              psixz = 0
              psiyy = 0
              psiyz = 0
              psizz = 0
           end where
        end if
        
        where ((x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= rad**2)
           gxx = 1
           gxy = 0
           gxz = 0
           gyy = 1
           gyz = 0
           gzz = 1
           kxx = 0
           kxy = 0
           kxz = 0
           kyy = 0
           kyz = 0
           kzz = 0
        end where
        
     end if
     
     if (overwrite_lapse /= 0) then
        
        where ((x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= rad**2)
           alp = 1
        end where
        
     end if
        
     if (overwrite_shift /= 0) then
        
        if (shift_state /= 0) then
           where ((x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= rad**2)
              betax = 0
              betay = 0
              betaz = 0
           end where
        end if
        
     end if
     
  end do
  
end subroutine NoExcision_Overwrite
