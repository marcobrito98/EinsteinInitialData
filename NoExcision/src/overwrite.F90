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
  
  CCTK_REAL, parameter :: zero=0
  CCTK_REAL :: dist2 (cctk_lsh(1), cctk_lsh(2), cctk_lsh(3))
  CCTK_REAL :: cx, cy, cz, rad, width
  integer   :: n
  
  do n = 1, num_regions
     
     cx = centre_x(n)
     cy = centre_y(n)
     cz = centre_z(n)
     rad = radius(n)
     width = smoothing_zone_width(n)
     
     dist2 = (x - cx)**2 + (y - cy)**2 + (z - cz)**2
     
     if (overwrite_geometry(n) /= 0) then
        
        if (conformal_state >= 1) then
           where (dist2 <= rad**2)
              psi = 1
           end where
        end if
        if (conformal_state >= 2) then
           where (dist2 <= rad**2)
              psix = 0
              psiy = 0
              psiz = 0
           end where
        end if
        if (conformal_state >= 3) then
           where (dist2 <= rad**2)
              psixx = 0
              psixy = 0
              psixz = 0
              psiyy = 0
              psiyz = 0
              psizz = 0
           end where
        end if
        
        where (dist2 <= rad**2)
           gxx = smooth (gxx,  Minkowski_scale(n), dist2)
           gxy = smooth (gxy,  zero              , dist2)
           gxz = smooth (gxz,  zero              , dist2)
           gyy = smooth (gyy,  Minkowski_scale(n), dist2)
           gyz = smooth (gyz,  zero              , dist2)
           gzz = smooth (gzz,  Minkowski_scale(n), dist2)
           kxx = smooth (kxx,  zero              , dist2)
           kxy = smooth (kxy,  zero              , dist2)
           kxz = smooth (kxz,  zero              , dist2)
           kyy = smooth (kyy,  zero              , dist2)
           kyz = smooth (kyz,  zero              , dist2)
           kzz = smooth (kzz,  zero              , dist2)
        end where
        
     end if
     
     if (overwrite_lapse(n) /= 0) then
        
        where (dist2 <= rad**2)
           alp = smooth (alp, lapse_scale(n), dist2)
        end where
        
     end if
        
     if (overwrite_shift(n) /= 0) then
        
        if (shift_state /= 0) then
           where (dist2 <= rad**2)
              betax = smooth (betax, zero, dist2)
              betay = smooth (betay, zero, dist2)
              betaz = smooth (betaz, zero, dist2)
           end where
        end if
        
     end if
     
  end do
  
contains
  
  elemental function smooth (oldval, goodval, dist2) result (newval)
    CCTK_REAL             :: newval
    CCTK_REAL, intent(in) :: oldval
    CCTK_REAL, intent(in) :: goodval
    CCTK_REAL, intent(in) :: dist2
    
    if (dist2 <= (rad-width)**2) then
       newval = goodval
    else if (width == 0 .or. dist2 >= rad**2) then
       newval = oldval
    else
       newval =   goodval &
            &   + (oldval - goodval) * (sqrt(dist2) - (rad-width)) / width
    end if
    
  end function smooth
  
end subroutine NoExcision_Overwrite
