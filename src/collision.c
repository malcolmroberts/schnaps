#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

void CollisionNumFlux(double wL[], double wR[], double *vnorm, double *flux)
{
  double vn =
    1 * vnorm[0] +
    0 * vnorm[1] +
    0 * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];
}

void CollisionBoundaryFlux(double x[3], double t, double wL[], double *vnorm,
			   double* flux)
{
  double wR[1];
  CollisionImposedData(x,t,wR);
  CollisionNumFlux(wL,wR,vnorm,flux);
}

void CollisionInitData(double x[3], double w[])
{
  double t=0;
  CollisionImposedData(x, t, w);
}

void CollisionImposedData(double x[3], double t, double w[])
{
  double vx =
    1 * x[0] +
    0 * x[1] +
    0 * x[2];
  double xx = vx - t;
  w[0] = cos(xx);
}
