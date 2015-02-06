#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>



void CollisionNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  //for(int i=0;i<_MV;i++){
  for(int i=0;i<1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    /* double vn = nel*_DV + */
    /*   _DV* gauss_lob_point[j+gauss_lob_offset[_DEG_V]]; */
    double vn = 1*vnorm[0];
    
    double vnp = vn>0 ? vn : 0;
    double vnm = vn-vnp;

    flux[i] = vnp * wL[i] + vnm * wR[i];
  }

};


void CollisionBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  CollisionImposedData(x,t,wR);
  CollisionNumFlux(wL,wR,vnorm,flux);
};


void CollisionInitData(double x[3],double w[]){

  double t=0;
  CollisionImposedData(x,t,w);

};



void CollisionImposedData(double x[3],double t,double w[]){

  double vx =
    1 * x[0] +
    0 * x[1] +
    0 * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};



