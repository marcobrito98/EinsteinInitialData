c param_defs.inc -- integer constants for decoded_exact_model
c $Header$

c
c For reasons explained in our param.ccl file, we decode the  exact_model
c parameter, and all our other string-valued parameters used in computing
c the stress-energy tensor, into integers.  This file contains #define
c definitions for those integers.
c

c
c For each parameter, the value 0 is deliberately *not* a legal value
c for the decoded integer, to help catch bugs where the decoded integer
c is not initialized properly.
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c ***** definitions for decoded_exact_model *****
c
c These *MUST* be distinct integers, and there must be precisely one
c definition for each value of  exact_model  in  param.ccl .  It is not
c necessary that they be in numerical order, that is just to make things
c look prettier.
c
c N.b. We are counting on the preprocessor being case-sensitive here,
c      since these same names with EXACT changed to Exact, are names
c      of subroutines for the individual metric types!
c

c Minkowski spacetime
#define EXACT__Minkowski		1
#define EXACT__Minkowski_shift		2
#define EXACT__Minkowski_funny		3
#define EXACT__Minkowski_gauge_wave	4
#define EXACT__Minkowski_conf_wave	5

c black hole spacetimes
#define EXACT__Schwarzschild_EF		10
#define EXACT__Schwarzschild_PG		11
#define EXACT__Schwarzschild_BL		12
#define EXACT__Schwarzschild_Novikov	13
#define EXACT__Kerr_BoyerLindquist	14
#define EXACT__Kerr_KerrSchild		15
#define EXACT__Schwarzschild_Lemaitre	16
#define EXACT__multi_BH			17
#define EXACT__Alvi			18
#define EXACT__Thorne_fakebinary	19

c cosmological spacetimes
#define EXACT__Lemaitre			50
CC this metric doesnt work, and has been moved to ./archive/
CC#define EXACT__Robertson_Walker	51
#define EXACT__de_Sitter		52
#define EXACT__de_Sitter_Lambda		53
#define EXACT__anti_de_Sitter_Lambda	54
#define EXACT__Bianchi_I		55
#define EXACT__Goedel			56
#define EXACT__Bertotti			57
#define EXACT__Kasner_like		58
#define EXACT__Kasner_axisymmetric	59
#define EXACT__Kasner_generalized	60
#define EXACT__Gowdy_wave		61
#define EXACT__Milne			62

c miscelaneous spacetimes
#define EXACT__boost_rotation_symmetric	80
#define EXACT__bowl			81
#define EXACT__constant_density_star	82

c math
#define EXACT__pi 3.14159265358979d0
