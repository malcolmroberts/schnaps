#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"



void CollisionNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  for(int i=0;i<_MV;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vn = (nel*_DV +
		 _DV* glop(_DEG_V,j))*vnorm[0];
    
    double vnp = vn>0 ? vn : 0;
    double vnm = vn-vnp;

    flux[i] = vnp * wL[i] + vnm * wR[i];
  }

};


void CollisionBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[_MV];
  CollisionImposedData(x,t,wR);
  CollisionNumFlux(wL,wR,vnorm,flux);
};


void CollisionInitData(double x[3],double w[]){

  double t=0;
  CollisionImposedData(x,t,w);

};



void CollisionImposedData(double x[3],double t,double w[]){

  for(int i=0;i<_MV;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vi = (nel*_DV +
		 _DV* glop(_DEG_V,j));
    w[i]=cos(x[0]-vi*t);
  }

};



