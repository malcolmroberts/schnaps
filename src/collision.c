#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_collision.h"



void Collision_Lagrangian_NumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  for(int i=0;i<_MV;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vn = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j))*vnorm[0];
    
    double vnp = vn>0 ? vn : 0;
    double vnm = vn-vnp;

    flux[i] = vnp * wL[i] + vnm * wR[i];
  }
  // do not change the potential !
  // and the electric field
  flux[_MV]=0;  // flux for phi
  flux[_MV+1]=0; // flux for E
  flux[_MV+2]=0; // flux for rho
  flux[_MV+3]=0; // flux for u
  flux[_MV+4]=0; // flux for p
  flux[_MV+5]=0; // flux for e ou T

};





void Collision_Lagrangian_BoundaryFlux(double x[3],double t,double wL[],double* vnorm,
				       double* flux){
  double wR[_MV+6];
  CollisionImposedData(x,t,wR);
  Collision_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};


void CollisionInitData(double x[3],double w[]){

  double t=0;
  CollisionImposedData(x,t,w);

};



void CollisionImposedData(double x[3],double t,double w[]){

  for(int i=0;i<_MV;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vi = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j));

    w[i]=Collision_ImposedKinetic_Data(x,t,vi);
  }
  // exact value of the potential
  // and electric field
  w[_MV]=x[0]*x[0];
  w[_MV+1]=1;
  w[_MV+2]=0; //rho init
  w[_MV+3]=0; // u init
  w[_MV+4]=0; // p init
  w[_MV+5]=0; // e ou T init

};




//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void CollisionSource(double* x,double t,double* w, double* source){

  double E=w[_MV+1]; // electric field
  double Md[_MV];
  double db[_MV];
  for(int iv=0;iv<_MV;iv++){
    Md[iv]=0;
    db[iv]=0;
  }
  for(int iv=0;iv<_MV+2;iv++){
    source[iv]=0;
  }
  // no source on the potential for the moment
  source[_MV]=0;
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<_DEG_V+1;kloc++){
      double omega=wglop(_DEG_V,kloc);
      int kpg=kloc+iel*_DEG_V;
      Md[kpg]+=omega*_DV;
      for(int iloc=0;iloc<_DEG_V+1;iloc++){
	int ipg=iloc+iel*_DEG_V;
	source[ipg]+=E*omega*w[kpg]*dlag(_DEG_V,iloc,kloc);
	if (iloc==kloc) db[ipg]+=E*omega*dlag(_DEG_V,iloc,kloc);
      }
    }
  }

  // upwinding
  if (E>0){
    source[_MV-1]-=E*w[_MV-1];
    db[_MV-1]-=E;
  }
  else {
    source[0]-=-E*w[0];
    db[0]-=-E;
  }

  for(int iv=0;iv<_MV;iv++){
    source[iv]/=Md[iv];
    //printf("%f ",db[iv]);
  }
  //printf("\n");
  //assert(1==2);
  

};


