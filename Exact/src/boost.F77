c This subroutine calculates the 4-metric and its inverse at an event,
c taking into account an optional Lorentz boost.
c $Header$
c
c The coordinates are
c   Cx(a) = Cactus $x^a$
c   Mx(a) = Model  $X^a$
c The 4-metrics are
c   Cgdd(a,b) = Cactus $g_{ab}$            Cguu(a,b) = Cactus $g^{ab}$
c   Mgdd(a,b) = Model  $g_{ab}$            Mguu(a,b) = Model  $g^{ab}$
c
c This file is copyright (c) 2003 by Jonathan Thornburg <jthorn@aei.mpg.de>.
c This file is covered by the GNU GPL license; see the files ../README
c and ../COPYING for details. 
c

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "param_defs.inc"

	subroutine Exact__metric(
     $     decoded_exact_model,
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz,
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $     gutt, gutx, guty, gutz,
     $     guxx, guyy, guzz, guxy, guyz, guxz,
     $     psi, rama)

	implicit none
	DECLARE_CCTK_FUNCTIONS
	DECLARE_CCTK_PARAMETERS

c input arguments
	CCTK_INT decoded_exact_model
	CCTK_REAL x, y, z, t

c output arguments
	CCTK_REAL gdtt, gdtx, gdty, gdtz,
     $		  gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $		  gutt, gutx, guty, gutz,
     $		  guxx, guyy, guzz, guxy, guyz, guxz,
     $            psi, rama

c intrinsic functions called
	CCTK_REAL sqrt

c static local variables describing Lorentz transformation
	logical   firstcall
        data      firstcall /.true./
	CCTK_REAL gamma
	CCTK_REAL vv(3), nn(3)
	CCTK_REAL parallel(3,3), perp(3,3)
	CCTK_REAL Cx_par(3), Cx_perp(3)
	CCTK_REAL partial_Mx_wrt_Cx(0:3,0:3)
	CCTK_REAL partial_Cx_wrt_Mx(0:3,0:3)
        save      firstcall
        save      gamma
        save      vv, nn
        save      parallel, perp
	save      Cx_par, Cx_perp
        save      partial_Mx_wrt_Cx
        save      partial_Cx_wrt_Mx

c coordinates and 4-metric
	CCTK_REAL Ct, Cx(3)
	CCTK_REAL Cgdd(0:3,0:3), Cguu(0:3,0:3)
	CCTK_REAL Mt, Mx(3)
	CCTK_REAL Mgdd(0:3,0:3), Mguu(0:3,0:3)

c misc temps
	CCTK_REAL vnorm, vnormsq
	CCTK_REAL delta_ij
	CCTK_REAL Cx_par_i, Cx_perp_i
	CCTK_REAL vdotCx
	CCTK_REAL Cgdd_ab, Cguu_ab
	character*100 warn_buffer

c flags, array indices, etc
	logical Tmunu_flag
	integer i, j
	integer Ca, Cb, MA, MB

c constants
	integer n
	parameter (n = 3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c optimized fast-path if no Lorentz boost
c
	if (       (boost_vx .eq. 0.0)
     $	     .and. (boost_vy .eq. 0.0)
     $	     .and. (boost_vz .eq. 0.0)       ) then
		call Exact__metric_for_model(
     $			decoded_exact_model,
     $			x, y, z, t,
     $			gdtt, gdtx, gdty, gdtz,
     $			gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $			gutt, gutx, guty, gutz,
     $			guxx, guyy, guzz, guxy, guyz, guxz,
     $			psi, Tmunu_flag,
     $			rama)
		return
	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c the rest of this function is the Lorentz-boost case:
c - Lorentz-transform Cactus coordinates --> Model coordinates
c - compute Model 4-metric and inverse at Model coordinates
c - tensor-transform 4-metric and inverse from Model coordinates
c   --> Cactus coordinates
c
c All the equations used are given in ../doc/documentation.tex
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c compute Lorentz transformation information on first call
c
	if (firstcall) then
		firstcall = .false.

c boost velocity
		vv(1) = boost_vx
		vv(2) = boost_vy
		vv(3) = boost_vz

c Lorentz gamma factor, unit vector in direction of boost velocity
		vnormsq = 0.0d0
			do 100 i = 1,n
			vnormsq = vnormsq + vv(i)*vv(i)
100			continue
		gamma = 1.0 / sqrt(1.0 - vnormsq)
		vnorm = sqrt(vnormsq)
			do 110 i = 1,n
			nn(i) = vv(i) / vnorm
110			continue

c projection operators parallel(*,*) and perp(*,*)
			do 210 j = 1,n
				do 200 i = 1,n
				parallel(i,j) = nn(i) * nn(j)
				if (i .eq. j) then
					delta_ij = 1.0d0
				   else
					delta_ij = 0.0d0
				endif
				perp(i,j) = delta_ij - parallel(i,j)
200				continue
210			continue

c partial derivatives of Model coordinates with respect to Cactus coordinates
		partial_Mx_wrt_Cx(0,0) = gamma
			do 300 i = 1,n
			partial_Mx_wrt_Cx(0,i) = -gamma*vv(i)
300			continue
			do 320 i = 1,n
			partial_Mx_wrt_Cx(i,0) = -gamma*vv(i)
				do 310 j=1,n
				partial_Mx_wrt_Cx(i,j)
     $					= gamma*parallel(i,j)
     $					  + perp(i,j)
310				continue
320			continue

c partial derivatives of Cactus coordinates with respect to Model coordinates
		partial_Cx_wrt_Mx(0,0) = gamma
			do 400 i = 1,n
			partial_Cx_wrt_Mx(0,i) = + gamma*vv(i)
400			continue
			do 420 i = 1,n
			partial_Cx_wrt_Mx(i,0) = + gamma*vv(i)
				do 410 j=1,n
				partial_Cx_wrt_Mx(i,j)
     $					= gamma*parallel(i,j)
     $					  + perp(i,j)
410				continue
420			continue

	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c compute flat-space components of Cx(*) parallel and perpendicular to vv(*)
c
	Ct    = t
	Cx(1) = x
	Cx(2) = y
	Cx(3) = z

			do 530 i=1,n
			Cx_par_i  = 0.0d0
			Cx_perp_i = 0.0d0
				do 520 j=1,n
				Cx_par_i  = Cx_par_i
     $					    + parallel(i,j)*Cx(j)
				Cx_perp_i = Cx_perp_i
     $					    + perp    (i,j)*Cx(j)
520				continue
			Cx_par (i) = Cx_par_i
			Cx_perp(i) = Cx_perp_i
530			continue

c
c Lorentz-transform the Cactus coordinate to get the Model coordinates
c
	vdotCx = 0.0
		do 600 i = 1,n
		vdotCx = vdotCx + vv(i)*Cx(i)
600		continue

	Mt = gamma * (Ct - vdotCx)
		do 610 i=1,n
		Mx(i) = gamma * (Cx_par(i) - vv(i)*Ct)
     $			+ Cx_perp(i)
610		continue

c
c compute the Model 4-metric and inverse 4-metric at the Model coordinates
c
	call Exact__metric_for_model(
     $			decoded_exact_model,
     $			Mx(1), Mx(2), Mx(3), Mt,
     $			Mgdd(0,0), Mgdd(0,1), Mgdd(0,2), Mgdd(0,3),
     $			Mgdd(1,1), Mgdd(2,2), Mgdd(3,3),
     $				Mgdd(1,2), Mgdd(2,3), Mgdd(1,3),
     $			Mguu(0,0), Mguu(0,1), Mguu(0,2), Mguu(0,3),
     $			Mguu(1,1), Mguu(2,2), Mguu(3,3),
     $				Mguu(1,2), Mguu(2,3), Mguu(1,3),
     $			psi, Tmunu_flag,
     $			rama)

	if (Tmunu_flag) then
		write (warn_buffer, '(a,i8,a,a)')
     $		      'exact_model = ', decoded_exact_model,
     $		      'sets the stress-energy tensor',
     $		      ' ==> we cannot Lorentz-boost it! :('
		call CCTK_WARN(0, warn_buffer)
	endif

c
c symmetrize the Model 4-metric and inverse 4-metric arrays
c (the Exact__metric_for_model() call only set the upper triangles)
c
	Mgdd(1,0) = Mgdd(0,1)
	Mgdd(2,0) = Mgdd(0,2)
	Mgdd(2,1) = Mgdd(1,2)
	Mgdd(3,0) = Mgdd(0,3)
	Mgdd(3,1) = Mgdd(1,3)
	Mgdd(3,2) = Mgdd(2,3)

	Mguu(1,0) = Mguu(0,1)
	Mguu(2,0) = Mguu(0,2)
	Mguu(2,1) = Mguu(1,2)
	Mguu(3,0) = Mguu(0,3)
	Mguu(3,1) = Mguu(1,3)
	Mguu(3,2) = Mguu(2,3)

c
c tensor-transorm (the upper triangle of) the 4-metric and inverse 4-metric
c
		do 730 Ca = 0,n
		do 720 Cb = Ca,n
		Cgdd_ab = 0.0d0
			do 710 Ma = 0,n
			do 700 Mb = 0,n
			Cgdd_ab = Cgdd_ab
     $				  +  Mgdd(Ma,Mb)
     $				     * partial_Mx_wrt_Cx(Ma,Ca)
     $				     * partial_Mx_wrt_Cx(Mb,Cb)
700			continue
710			continue
		Cgdd(Ca,Cb) = Cgdd_ab
720		continue
730		continue

		do 830 Ca = 0,n
		do 820 Cb = Ca,n
		Cguu_ab = 0.0d0
			do 810 Ma = 0,n
			do 800 Mb = 0,n
			Cguu_ab = Cguu_ab
     $				  +  Mguu(Ma,Mb)
     $				     * partial_Cx_wrt_Mx(Ca,Ma)
     $				     * partial_Cx_wrt_Mx(Cb,Mb)
800			continue
810			continue
		Cguu(Ca,Cb) = Cguu_ab
820		continue
830		continue

c
c unpack the Cactus-coordinates 4-metric and inverse 4-metric
c into the corresponding output arguments
c
	gdtt = Cgdd(0,0)
	gdtx = Cgdd(0,1)
	gdty = Cgdd(0,2)
	gdtz = Cgdd(0,3)
	gdxx = Cgdd(1,1)
	gdxy = Cgdd(1,2)
	gdxz = Cgdd(1,3)
	gdyy = Cgdd(2,2)
	gdyz = Cgdd(2,3)
	gdzz = Cgdd(3,3)

	gutt = Cguu(0,0)
	gutx = Cguu(0,1)
	guty = Cguu(0,2)
	gutz = Cguu(0,3)
	guxx = Cguu(1,1)
	guxy = Cguu(1,2)
	guxz = Cguu(1,3)
	guyy = Cguu(2,2)
	guyz = Cguu(2,3)
	guzz = Cguu(3,3)

	return
	end
