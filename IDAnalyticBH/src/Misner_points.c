 /*@@
   @file      Misner_points.c
   @date      
   @author    Steve Brandt
   @desc 
      	      This calculates the conformal factor for nbh black holes, 
              with naked mass m0 = 2 csch(mu) each, and placed on a circle in the
              xy plane around the origin, of radius coth(mu). 
              One of them sits on the positive x axis, the others are evenly spaced.
              Naked mass here corresponds to the term m0 / (2 |r - r0|) in the expansion.
   @enddesc 
 @@*/

#include "cctk.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

void Misner_init(int n, CCTK_REAL mu, int terms);
void MisnerEvalPsi(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *res);

static int nbholes;

/* Basic data about a brill-lindquist black hole term. */
struct bhole {
  CCTK_REAL x,y;
  CCTK_REAL mass;

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

static CCTK_REAL csch(CCTK_REAL theta) {
  return 1.0/sinh(theta); 
}
static CCTK_REAL coth(CCTK_REAL theta) {
  return cosh(theta)/sinh(theta); 
}

 /*@@
   @routine    fill_iso
   @date       
   @author     Steve Brandt
   @desc 
     	       Isometrize black hole a1 through hole a2
   @enddesc 
   @calls      fill_iso   
   @history 
 
   @endhistory 

@@*/

static void iso(struct bhole *a1, struct bhole *a2, struct bhole *a3) 
{
  CCTK_REAL rad,radtwo;
  radtwo=(
      (a1->x - a2->x)*(a1->x - a2->x)+
      (a1->y - a2->y)*(a1->y - a2->y)
    );
  rad=sqrt(radtwo);
  a3->mass = a1->mass*a2->mass/rad;
  a3->x = a2->x+(a2->mass*a2->mass)*(a1->x - a2->x)/radtwo;
  a3->y = a2->y+(a2->mass*a2->mass)*(a1->y - a2->y)/radtwo;
}

 /*@@
   @routine    fill_iso
   @date       
   @author     Steve Brandt
   @desc 
     	       Fills in the iso structure of a given black hole.  
               Applies recursively to the number of terms desired.
   @enddesc 
   @calls      fill_iso   
   @history 
 
   @endhistory 

@@*/

static void fill_iso(struct bhole *b, int n) 
{
  int i,j;
  if(n==0) 
  {
    b->isos = 0;
    return;
  }
  b->isos = (struct bhole *)malloc(sizeof(struct bhole)*(nbholes-1));
  assert(b->isos != 0);
  for(j=0, i=0;i<nbholes;i++) 
  {
    if(i != b->i) {
      iso(b,&bholes[i],&b->isos[j]);
      b->isos[j].i = i;
      fill_iso(&b->isos[j],n-1);
      j++;
    }
  }
}

 /*@@
   @routine    Misner_init
   @date       
   @author     Steve Brandt
   @desc 
      Initialises the black holes then makes the isometry black holes
   @enddesc 
   @calls     
   @history 
 
   @endhistory 

@@*/

void Misner_init(int n, CCTK_REAL mu, int terms) 
{

  int i;
  CCTK_REAL pi,ang;

  assert((nbholes=n) < MAXBHOLES);

  pi = 4.0*atan(1.);

  ang = 2.*pi/(n);

  for(i=0;i<n;i++) 
  {
    bholes[i].x = coth(mu)*cos(ang*i);
    bholes[i].y = coth(mu)*sin(ang*i);
    bholes[i].mass = csch(mu);
    bholes[i].i = i;
    bholes[i].isos = 0;
  }

  for(i=0;i<n;i++)
    fill_iso(&bholes[i],terms);

}


 /*@@
   @routine    eval_bh_psi
   @date       
   @author     Steve Brandt
   @desc 
      
   @enddesc 
   @calls      eval_bh_psi   
   @history 
 
   @endhistory 

@@*/

static CCTK_REAL eval_bh_psi(struct bhole *b, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) 
{
  int i;
  CCTK_REAL res;
  res = 0.0;
  if(b->isos != 0) 
  {
    for(i=0;i<nbholes-1;i++) 
    {
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


 /*@@
   @routine    MisnerEvalPsi
   @date       
   @author     Steve Brandt
   @desc 
      Evaluate psi at a point
   @enddesc 
   @calls      eval_bh_psi
   @history 
 
   @endhistory 

@@*/

void MisnerEvalPsi(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *res) 
{
  int i;
  *res = 1;
  for(i=0;i<nbholes;i++)
    *res += eval_bh_psi(&bholes[i],x,y,z);
}
