C Here we added the matter variables for several of the metrics
C you can find in "src" directory through the components of the
C stress-eergy tensor. Being different in specific cases it is
C necessary to check first which metric is running.  
C Author : Dumitru Vulcanov (Timisoara, Romania)


C Kasner like metric

	if (CCTK_Equals(exactmodel,"kasner").eq.1) then

	Ttt = Ttt - kasner_q*(2-3*kasner_q)/(CCTK_time**2)
	
C DeSitter metric

	elseif (CCTK_Equals(exactmodel,"desitter").eq.1) then

	Ttt =  Ttt - desitt_b/(CCTK_time**2)

C Robertson Walker cosmology

	elseif (CCTK_Equals(exactmodel,"rob-wal").eq.1) then

          rr2  = x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2
          aha  = desitt_b/(3*desitt_a**2*CCTK_time**2)

	Ttt =  Ttt + desitt_b/(CCTK_time**2)/(desitt_a**4)

        Txx =  Txx + 1.d0+ (aha - 1.d0)*x(i,j,k)**2/rr2

        Tyy =  Tyy + 1.d0+ (aha - 1.d0)*y(i,j,k)**2/rr2

        Tzz =  Tzz + 1.d0+ (aha - 1.d0)*z(i,j,k)**2/rr2

        Txy =  Txy + (aha-1.d0)*x(i,j,k)*y(i,j,k)/rr2

        Txz =  Txz + (aha-1.d0)*x(i,j,k)*z(i,j,k)/rr2

        Tyz =  Tyz + (aha-1.d0)*y(i,j,k)*z(i,j,k)/rr2

C Schwarzschild static star 

	elseif (CCTK_Equals(exactmodel,"starschwarz").eq.1) then

          rr  = dsqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)

	if (rr.le.starschwarz_r) then


	  sum = gxx(i,j,k)*betax(i,j,k)**2+gxx(i,j,k)*betay(i,j,k)**2 +
     &       gzz(i,j,k)*betaz(i,j,k)**2 +
     &       2*gxy(i,j,k)*betax(i,j,k)*betay(i,j,k) +  
     &       2*gxz(i,j,k)*betax(i,j,k)*betaz(i,j,k) +
     &       2*gyz(i,j,k)*betay(i,j,k)*betaz(i,j,k) 
     &       -alp(i,j,k)**2


   	  unu = 3*starschwarz_m/(4*pi*starschwarz_m**3)

	  doi = dsqrt(1-2*starschwarz_m/starschwarz_r)

	  trei= dsqrt(1-2*starschwarz_m*rr**2/(starschwarz_r**3))
	
	Ttt = Ttt +unu*alp(i,j,k)**2 - (unu*(3*doi/2-trei/2)**2)*sum

	Txx = Txx - (unu*(doi-trei)/(trei-3*unu))*gxx(i,j,k) 

	Tyy = Tyy - (unu*(doi-trei)/(trei-3*unu))*gyy(i,j,k) 

	Tzz = Tzz - (unu*(doi-trei)/(trei-3*unu))*gzz(i,j,k) 


	else
	Ttt = Ttt
	Txx = Txx
	Tyy = Tyy
	Tzz = Tzz

	end if 


        end if 