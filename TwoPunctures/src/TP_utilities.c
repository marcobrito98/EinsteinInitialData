/* TwoPunctures:  File  "utilities.c"*/

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "TP_utilities.h"

#include "cctk_Functions.h"

/*---------------------------------------------------------------------------*/
int *
ivector (long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v = (int *) malloc ((size_t) ((nh - nl + 1 + NR_END) * sizeof (int)));
  if (!v)
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure in ivector()");
  return v - nl + NR_END;
}

/*---------------------------------------------------------------------------*/
CCTK_REAL *
dvector (long nl, long nh)
/* allocate a CCTK_REAL vector with subscript range v[nl..nh] */
{
  CCTK_REAL *v;

  v = (CCTK_REAL *) malloc ((size_t) ((nh - nl + 1 + NR_END) * sizeof (CCTK_REAL)));
  if (!v)
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure in dvector()");
  return v - nl + NR_END;
}

/*---------------------------------------------------------------------------*/
int **
imatrix (long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc ((size_t) ((nrow + NR_END) * sizeof (int *)));
  if (!m)
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl] = (int *) malloc ((size_t) ((nrow * ncol + NR_END) * sizeof (int)));
  if (!m[nrl])
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/*---------------------------------------------------------------------------*/
CCTK_REAL **
dmatrix (long nrl, long nrh, long ncl, long nch)
/* allocate a CCTK_REAL matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  CCTK_REAL **m;

  /* allocate pointers to rows */
  m = (CCTK_REAL **) malloc ((size_t) ((nrow + NR_END) * sizeof (CCTK_REAL *)));
  if (!m)
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] =
    (CCTK_REAL *) malloc ((size_t) ((nrow * ncol + NR_END) * sizeof (CCTK_REAL)));
  if (!m[nrl])
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/*---------------------------------------------------------------------------*/
CCTK_REAL ***
d3tensor (long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a CCTK_REAL 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  CCTK_REAL ***t;

  /* allocate pointers to pointers to rows */
  t = (CCTK_REAL ***) malloc ((size_t) ((nrow + NR_END) * sizeof (CCTK_REAL **)));
  if (!t)
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] =
    (CCTK_REAL **)
    malloc ((size_t) ((nrow * ncol + NR_END) * sizeof (CCTK_REAL *)));
  if (!t[nrl])
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] =
    (CCTK_REAL *)
    malloc ((size_t) ((nrow * ncol * ndep + NR_END) * sizeof (CCTK_REAL)));
  if (!t[nrl][ncl])
    CCTK_WARN (CCTK_WARN_ABORT, "allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
  {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

/*--------------------------------------------------------------------------*/
void
free_ivector (int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free ((FREE_ARG) (v + nl - NR_END));
}

/*--------------------------------------------------------------------------*/
void
free_dvector (CCTK_REAL *v, long nl, long nh)
/* free a CCTK_REAL vector allocated with dvector() */
{
  free ((FREE_ARG) (v + nl - NR_END));
}

/*--------------------------------------------------------------------------*/
void
free_imatrix (int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free ((FREE_ARG) (m[nrl] + ncl - NR_END));
  free ((FREE_ARG) (m + nrl - NR_END));
}

/*--------------------------------------------------------------------------*/
void
free_dmatrix (CCTK_REAL **m, long nrl, long nrh, long ncl, long nch)
/* free a CCTK_REAL matrix allocated by dmatrix() */
{
  free ((FREE_ARG) (m[nrl] + ncl - NR_END));
  free ((FREE_ARG) (m + nrl - NR_END));
}

/*--------------------------------------------------------------------------*/
void
free_d3tensor (CCTK_REAL ***t, long nrl, long nrh, long ncl, long nch,
	       long ndl, long ndh)
/* free a CCTK_REAL f3tensor allocated by f3tensor() */
{
  free ((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
  free ((FREE_ARG) (t[nrl] + ncl - NR_END));
  free ((FREE_ARG) (t + nrl - NR_END));
}

/*--------------------------------------------------------------------------*/
int
minimum2 (int i, int j)
{
  int result = i;
  if (j < result)
    result = j;
  return result;
}

/*-------------------------------------------------------------------------*/
int
minimum3 (int i, int j, int k)
{
  int result = i;
  if (j < result)
    result = j;
  if (k < result)
    result = k;
  return result;
}

/*--------------------------------------------------------------------------*/
int
maximum2 (int i, int j)
{
  int result = i;
  if (j > result)
    result = j;
  return result;
}

/*--------------------------------------------------------------------------*/
int
maximum3 (int i, int j, int k)
{
  int result = i;
  if (j > result)
    result = j;
  if (k > result)
    result = k;
  return result;
}

/*--------------------------------------------------------------------------*/
int
pow_int (int mantisse, int exponent)
{
  int i, result = 1;

  for (i = 1; i <= exponent; i++)
    result *= mantisse;

  return result;
}

/*--------------------------------------------------------------------------*/
#if 0
CCTK_REAL
atanh (CCTK_REAL x)
{
  return 0.5 * log ((1 + x) / (1 - x));
}

/*--------------------------------------------------------------------------*/
CCTK_REAL
asinh (CCTK_REAL x)
{
  return log (x + sqrt (1 + x * x));
}

/*--------------------------------------------------------------------------*/
CCTK_REAL
acosh (CCTK_REAL x)
{
  return log (x + sqrt (x * x - 1));
}
#endif

/*--------------------------------------------------------------------------*/
dcomplex
Cadd (dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;
  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Csub (dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;
  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Cmul (dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r * b.r - a.i * b.i;
  c.i = a.i * b.r + a.r * b.i;
  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
RCmul (CCTK_REAL x, dcomplex a)
{
  dcomplex c;
  c.r = x * a.r;
  c.i = x * a.i;
  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Cdiv (dcomplex a, dcomplex b)
{
  dcomplex c;
  CCTK_REAL r, den;
  if (fabs (b.r) >= fabs (b.i))
  {
    r = b.i / b.r;
    den = b.r + r * b.i;
    c.r = (a.r + r * a.i) / den;
    c.i = (a.i - r * a.r) / den;
  }
  else
  {
    r = b.r / b.i;
    den = b.i + r * b.r;
    c.r = (a.r * r + a.i) / den;
    c.i = (a.i * r - a.r) / den;
  }
  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Complex (CCTK_REAL re, CCTK_REAL im)
{
  dcomplex c;
  c.r = re;
  c.i = im;
  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Conjg (dcomplex z)
{
  dcomplex c;
  c.r = z.r;
  c.i = -z.i;
  return c;
}

/*--------------------------------------------------------------------------*/
CCTK_REAL
Cabs (dcomplex z)
{
  CCTK_REAL x, y, ans, temp;
  x = fabs (z.r);
  y = fabs (z.i);
  if (x == 0.0)
    ans = y;
  else if (y == 0.0)
    ans = x;
  else if (x > y)
  {
    temp = y / x;
    ans = x * sqrt (1.0 + temp * temp);
  }
  else
  {
    temp = x / y;
    ans = y * sqrt (1.0 + temp * temp);
  }
  return ans;
}

/*--------------------------------------------------------------------------*/
dcomplex
Csqrt (dcomplex z)
{
  dcomplex c;
  CCTK_REAL x, y, w, r;
  if ((z.r == 0.0) && (z.i == 0.0))
  {
    c.r = 0.0;
    c.i = 0.0;
    return c;
  }
  else
  {
    x = fabs (z.r);
    y = fabs (z.i);
    if (x >= y)
    {
      r = y / x;
      w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + r * r)));
    }
    else
    {
      r = x / y;
      w = sqrt (y) * sqrt (0.5 * (r + sqrt (1.0 + r * r)));
    }
    if (z.r >= 0.0)
    {
      c.r = w;
      c.i = z.i / (2.0 * w);
    }
    else
    {
      c.i = (z.i >= 0) ? w : -w;
      c.r = z.i / (2.0 * c.i);
    }
    return c;
  }
}

/*--------------------------------------------------------------------------*/
dcomplex
Cexp (dcomplex z)
{
  dcomplex c;
  CCTK_REAL exp_r = exp (z.r);

  c.r = exp_r * cos (z.i);
  c.i = exp_r * sin (z.i);

  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Clog (dcomplex z)
{
  dcomplex c;

  c.r = 0.5 * log (z.r * z.r + z.i * z.i);
  c.i = atan2 (z.i, z.r);

  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Csin (dcomplex z)
{
  dcomplex c;

  c.r = sin (z.r) * cosh (z.i);
  c.i = cos (z.r) * sinh (z.i);

  return c;
}
/*--------------------------------------------------------------------------*/

dcomplex
Ccos (dcomplex z)
{
  dcomplex c;

  c.r = cos (z.r) * cosh (z.i);
  c.i = -sin (z.r) * sinh (z.i);

  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Ctan (dcomplex z)
{
  return Cdiv (Csin (z), Ccos (z));
}

/*--------------------------------------------------------------------------*/
dcomplex
Ccot (dcomplex z)
{
  return Cdiv (Ccos (z), Csin (z));
}

/*--------------------------------------------------------------------------*/
dcomplex
Csinh (dcomplex z)
{
  dcomplex c;

  c.r = sinh (z.r) * cos (z.i);
  c.i = cosh (z.r) * sin (z.i);

  return c;
}
/*--------------------------------------------------------------------------*/

dcomplex
Ccosh (dcomplex z)
{
  dcomplex c;

  c.r = cosh (z.r) * cos (z.i);
  c.i = sinh (z.r) * sin (z.i);

  return c;
}

/*--------------------------------------------------------------------------*/
dcomplex
Ctanh (dcomplex z)
{
  return Cdiv (Csinh (z), Ccosh (z));
}

/*--------------------------------------------------------------------------*/
dcomplex
Ccoth (dcomplex z)
{
  return Cdiv (Ccosh (z), Csinh (z));
}

/*--------------------------------------------------------------------------*/
void
chebft_Zeros (CCTK_REAL u[], int n, int inv)
    /* eq. 5.8.7 and 5.8.8 at x = (5.8.4) of 2nd edition C++ NR */
{
  int k, j, isignum;
  CCTK_REAL fac, sum, Pion, *c;

  c = dvector (0, n);
  Pion = Pi / n;
  if (inv == 0)
  {
    fac = 2.0 / n;
    isignum = 1;
    for (j = 0; j < n; j++)
    {
      sum = 0.0;
      for (k = 0; k < n; k++)
	sum += u[k] * cos (Pion * j * (k + 0.5));
      c[j] = fac * sum * isignum;
      isignum = -isignum;
    }
  }
  else
  {
    for (j = 0; j < n; j++)
    {
      sum = -0.5 * u[0];
      isignum = 1;
      for (k = 0; k < n; k++)
      {
	sum += u[k] * cos (Pion * (j + 0.5) * k) * isignum;
	isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
#if 0
    if (fabs(c[j]) < 5.e-16)
      u[j] = 0.0;
    else
#endif
      u[j] = c[j];
  free_dvector (c, 0, n);
}

/* --------------------------------------------------------------------------*/

void
chebft_Extremes (CCTK_REAL u[], int n, int inv)
    /* eq. 5.8.7 and 5.8.8 at x = (5.8.5) of 2nd edition C++ NR */
{
  int k, j, isignum, N = n - 1;
  CCTK_REAL fac, sum, PioN, *c;

  c = dvector (0, N);
  PioN = Pi / N;
  if (inv == 0)
  {
    fac = 2.0 / N;
    isignum = 1;
    for (j = 0; j < n; j++)
    {
      sum = 0.5 * (u[0] + u[N] * isignum);
      for (k = 1; k < N; k++)
	sum += u[k] * cos (PioN * j * k);
      c[j] = fac * sum * isignum;
      isignum = -isignum;
    }
    c[N] = 0.5 * c[N];
  }
  else
  {
    for (j = 0; j < n; j++)
    {
      sum = -0.5 * u[0];
      isignum = 1;
      for (k = 0; k < n; k++)
      {
	sum += u[k] * cos (PioN * j * k) * isignum;
	isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
    u[j] = c[j];
  free_dvector (c, 0, N);
}

/* --------------------------------------------------------------------------*/

void
chder (CCTK_REAL *c, CCTK_REAL *cder, int n)
{
  int j;

  cder[n] = 0.0;
  cder[n - 1] = 0.0;
  for (j = n - 2; j >= 0; j--)
    cder[j] = cder[j + 2] + 2 * (j + 1) * c[j + 1];
}

/* --------------------------------------------------------------------------*/
CCTK_REAL
chebev (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c[], int m, CCTK_REAL x)
    /* eq. 5.8.11 of C++ NR (2nd ed) */
{
  int j;
  CCTK_REAL djp2, djp1, dj; /* d_{j+2}, d_{j+1} and d_j */
  CCTK_REAL y;

  /* rescale input to lie within [-1,1] */
  y = 2*(x - 0.5*(b+a))/(b-a);

  dj = djp1 = 0;
  for(j = m-1 ; j >= 1; j--)
  { 
    /* advance the coefficients */
    djp2 = djp1; 
    djp1 = dj;
    dj   = 2*y*djp1 - djp2 + c[j];
  }

  return y*dj - djp1 + 0.5*c[0];
}

/* --------------------------------------------------------------------------*/
void
fourft (CCTK_REAL *u, int N, int inv)
    /* a (slow) Fourier transform, seems to be just eq. 12.1.6 and 12.1.9 of C++ NR (2nd ed) */
{
  int l, k, iy, M;
  CCTK_REAL x, x1, fac, Pi_fac, *a, *b;

  M = N / 2;
  a = dvector (0, M);
  b = dvector (1, M);		/* Actually: b=vector(1,M-1) but this is problematic if M=1*/
  fac = 1. / M;
  Pi_fac = Pi * fac;
  if (inv == 0)
  {
    for (l = 0; l <= M; l++)
    {
      a[l] = 0;
      if (l > 0 && l < M)
	b[l] = 0;
      x1 = Pi_fac * l;
      for (k = 0; k < N; k++)
      {
	x = x1 * k;
	a[l] += fac * u[k] * cos (x);
	if (l > 0 && l < M)
	  b[l] += fac * u[k] * sin (x);
      }
    }
    u[0] = a[0];
    u[M] = a[M];
    for (l = 1; l < M; l++)
    {
      u[l] = a[l];
      u[l + M] = b[l];
    }
  }
  else
  {
    a[0] = u[0];
    a[M] = u[M];
    for (l = 1; l < M; l++)
    {
      a[l] = u[l];
      b[l] = u[M + l];
    }
    iy = 1;
    for (k = 0; k < N; k++)
    {
      u[k] = 0.5 * (a[0] + a[M] * iy);
      x1 = Pi_fac * k;
      for (l = 1; l < M; l++)
      {
	x = x1 * l;
	u[k] += a[l] * cos (x) + b[l] * sin (x);
      }
      iy = -iy;
    }
  }
  free_dvector (a, 0, M);
  free_dvector (b, 1, M);
}

/* -----------------------------------------*/
void
fourder (CCTK_REAL u[], CCTK_REAL du[], int N)
{
  int l, M, lpM;

  M = N / 2;
  du[0] = 0.;
  du[M] = 0.;
  for (l = 1; l < M; l++)
  {
    lpM = l + M;
    du[l] = u[lpM] * l;
    du[lpM] = -u[l] * l;
  }
}

/* -----------------------------------------*/
void
fourder2 (CCTK_REAL u[], CCTK_REAL d2u[], int N)
{
  int l, l2, M, lpM;

  d2u[0] = 0.;
  M = N / 2;
  for (l = 1; l <= M; l++)
  {
    l2 = l * l;
    lpM = l + M;
    d2u[l] = -u[l] * l2;
    if (l < M)
      d2u[lpM] = -u[lpM] * l2;
  }
}

/* ----------------------------------------- */
CCTK_REAL
fourev (CCTK_REAL *u, int N, CCTK_REAL x)
{
  int l, M = N / 2;
  CCTK_REAL xl, result;

  result = 0.5 * (u[0] + u[M] * cos (x * M));
  for (l = 1; l < M; l++)
  {
    xl = x * l;
    result += u[l] * cos (xl) + u[M + l] * sin (xl);
  }
  return result;
}

/* ------------------------------------------------------------------------*/
CCTK_REAL
norm1 (CCTK_REAL *v, int n)
{
  int i;
  CCTK_REAL result = -1;

  for (i = 0; i < n; i++)
    if (fabs (v[i]) > result)
      result = fabs (v[i]);

  return result;
}

/* -------------------------------------------------------------------------*/
CCTK_REAL
norm2 (CCTK_REAL *v, int n)
{
  int i;
  CCTK_REAL result = 0;

  for (i = 0; i < n; i++)
    result += v[i] * v[i];

  return sqrt (result);
}

/* -------------------------------------------------------------------------*/
CCTK_REAL
scalarproduct (CCTK_REAL *v, CCTK_REAL *w, int n)
{
  int i;
  CCTK_REAL result = 0;

  for (i = 0; i < n; i++)
    result += v[i] * w[i];

  return result;
}

/* -------------------------------------------------------------------------*/
