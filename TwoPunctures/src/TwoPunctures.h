// TwoPunctures:  File  "TwoPunctures.h"

#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

typedef struct DERIVS
{
  double *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

/*
Files of "TwoPunctures":
	TwoPunctures.c
	FuncAndJacobian.c
	CoordTransf.c
	Equations.c
	Newton.c
	utilities.c (see utilities.h)
**************************
*/

// Routines in  "TwoPunctures.c"
double TestSolution (double A, double B, double X, double R, double phi);
void TestVector_w (double *par, int nvar, int n1, int n2, int n3, double *w);

// Routines in  "FuncAndJacobian.c"
int Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void allocate_derivs (derivs * v, int n);
void free_derivs (derivs * v, int n);
void Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs v);
void F_of_v (int nvar, int n1, int n2, int n3, derivs v,
	     double *F, derivs u);
void J_times_dv (int nvar, int n1, int n2, int n3, derivs dv,
		 double *Jdv, derivs u);
void JFD_times_dv (int i, int j, int k, int nvar, int n1,
		   int n2, int n3, derivs dv, derivs u, double *values);
void SetMatrix_JFD (int nvar, int n1, int n2, int n3,
		    derivs u, int *ncols, int **cols, double **Matrix);
double PunctEvalAtArbitPosition (double *v, double A, double B, double phi,
				 int n1, int n2, int n3);
void calculate_derivs (int i, int j, int k, int ivar, int nvar, int n1,
		       int n2, int n3, derivs v, derivs vv);
double interpol (double a, double b, double c, derivs v);
double PunctIntPolAtArbitPosition (int ivar, int nvar, int n1,
				   int n2, int n3, derivs v, double x,
				   double y, double z);

// Routines in  "CoordTransf.c"
void AB_To_XR (int nvar, double A, double B, double *X,
	       double *R, derivs U);
void C_To_c (int nvar, double X, double R, double *x,
	     double *r, derivs U);
void rx3_To_xyz (int nvar, double x, double r, double phi, double *y,
		 double *z, derivs U);

// Routines in  "Equations.c"
double BY_KKofxyz (double x, double y, double z);
void BY_Aijofxyz (double x, double y, double z, double Aij[3][3]);
void NonLinEquations (double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs U, double *values);
void LinEquations (double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs dU, derivs U, double *values);

// Routines in  "Newton.c"
void resid (double *res, int ntotal, double *dv, double *rhs,
	    int *ncols, int **cols, double **JFD);
void LineRelax_al (double *dv, int j, int k, int nvar, int n1, int n2, int n3,
		   double *rhs, int *ncols, int **cols, double **JFD);
void LineRelax_be (double *dv, int i, int k, int nvar, int n1, int n2, int n3,
		   double *rhs, int *ncols, int **cols, double **JFD);
void relax (double *dv, int nvar, int n1, int n2, int n3, double *rhs,
	    int *ncols, int **cols, double **JFD);
void TestRelax (int nvar, int n1, int n2, int n3, derivs v,
		double *dv);
int bicgstab (int nvar, int n1, int n2, int n3, derivs v,
	      derivs dv, int output, int itmax, double tol, double *normres);
void Newton (int nvar, int n1, int n2, int n3, derivs v,
	     double tol, int itmax);


/* 
 27: -1.325691774825335e-03
 37: -1.325691778944117e-03
 47: -1.325691778942711e-03
 
 17: -1.510625972641537e-03
 21: -1.511443006977708e-03
 27: -1.511440785153687e-03
 37: -1.511440809549005e-03
 39: -1.511440809597588e-03
 */
