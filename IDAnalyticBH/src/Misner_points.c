/* This code written by Steve Brandt. */

/* This calculates the conformal factor for nbh black holes, 
with naked mass m0 = 2 csch(mu) each, and placed on a circle in the
xy plane around the origin, of radius coth(mu). 
One of them sits on the positive x axis, the others are evenly spaced.
Naked mass here corresponds to the term m0 / (2 |r - r0|) in the expansion.
For nbh =2, these are the same masses as in the routine misner.F, 
with the difference that here two BHs sit on the positive and negative
x axis, there on the z-axis. */

#include "cctk.h"

#include <stdio.h>
#include <assert.h>
#include <math.h>

int nbholes;
double mu;

/* Basic data about a brill-lindquist black hole term. */
struct bhole {
  double x,y;
  double mass;

  /* i gives either the number of the seed
     black hole we are starting with, or
     the number of the seed black hole that
     was used to isometrize this term. */
  int i;
  struct bhole *isos;
};

#define MAXBHOLES 10
/* The seed black holes. */
struct bhole bholes[MAXBHOLES];

double csch(double mu) {
  return 1.0/sinh(mu); 
}
double coth(double mu) {
  return cosh(mu)/sinh(mu); 
}

/* Isometrize black hole a1 through hole a2 */
void iso(struct bhole *a1, struct bhole *a2, struct bhole *a3) {
  double rad,radtwo;
  radtwo=(
      (a1->x - a2->x)*(a1->x - a2->x)+
      (a1->y - a2->y)*(a1->y - a2->y)
    );
  rad=sqrt(radtwo);
  a3->mass = a1->mass*a2->mass/rad;
  a3->x = a2->x+(a2->mass*a2->mass)*(a1->x - a2->x)/radtwo;
  a3->y = a2->y+(a2->mass*a2->mass)*(a1->y - a2->y)/radtwo;
}

/* Fills in the iso structure of a given black hole.  Applies
   recursively to the number of terms desired. */
void fill_iso(struct bhole *b,int n) {
  int i,j;
  if(n==0) {
    b->isos = 0;
    return;
  }
  b->isos = (struct bhole *)malloc(sizeof(struct bhole)*(nbholes-1));
  assert(b->isos != 0);
  for(j=0, i=0;i<nbholes;i++) {
    if(i != b->i) {
      iso(b,&bholes[i],&b->isos[j]);
      b->isos[j].i = i;
      fill_iso(&b->isos[j],n-1);
      j++;
    }
  }
}

/* Initializes the black holes, then makes the
   isometry black holes. */
void FORTRAN_NAME(nmisner_init)(int *n, double *mu,int *terms) 
{
  int i;
  double pi,ang;
  assert((nbholes=*n) < MAXBHOLES);
  pi = 4.0*atan(1.);
  ang = 2.*pi/(*n);
  for(i=0;i<*n;i++) {
    bholes[i].x = coth(*mu)*cos(ang*i);
    bholes[i].y = coth(*mu)*sin(ang*i);
    bholes[i].mass = csch(*mu);
    bholes[i].i = i;
    bholes[i].isos = 0;
  }
  for(i=0;i<*n;i++)
    fill_iso(&bholes[i],*terms);
}

double eval_bh_psi(struct bhole *b,double x,double y,double z) {
  int i;
  double res;
  res = 0.0;
  if(b->isos != 0) {
    for(i=0;i<nbholes-1;i++) {
      res += eval_bh_psi(&b->isos[i],x,y,z);
    }
  }
  res += b->mass/sqrt(
      (x - b->x)*(x - b->x)+
      (y - b->y)*(y - b->y)+
      z*z
    );
  return res;
}

/* use this function to evaluate psi at a point. */
void FORTRAN_NAME(nmisner_eval_psi)(double *x,double *y,double *z,double *res) {
  int i;
  *res = 1;
  for(i=0;i<nbholes;i++)
    *res += eval_bh_psi(&bholes[i],*x,*y,*z);
}
