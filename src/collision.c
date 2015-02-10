#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"


void Collision_Lagrangian_NumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
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


void Collision_Lagrangian_BoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[_MV];
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

    double vi = (nel*_DV +
		 _DV* glop(_DEG_V,j));
    w[i]=cos(x[0]-vi*t);
  }

};


double Collision_ImposedKinetic_Data(double x[3],double t,double v){
  double f;
  f=cos(x[0]-v*t);
  return f;
};

double L2_Kinetic_error(Field* f){

  double error=0;
  double error_space=0;
  double moy=0; // mean value
  double moy_space=0;

  for(int i=0;i<_MV;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vn = (nel*_DV +
		 _DV* glop(_DEG_V,j));
    
    double weight = wglop(_DEG_V,j);

    for (int ie=0;ie<f->macromesh.nbelems;ie++){
      // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],xphy[3],wpg;
      double dtau[3][3],codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      double w[f->model.m];
      double wex=0;
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      wex=Collision_ImposedKinetic_Data(xphy,f->tnow,vn);
      for(int iv=0;iv<f->model.m;iv++){
        error_space+=pow(w[iv]-wex,2)*wpg*det;
        moy_space+=pow(w[iv],2)*wpg*det;
      }
    }
  }
    error=error+weight*error_space;
    moy=moy+weight*moy_space;

    error_space=0;
    moy_space=0;
    
  }
  return sqrt(error)/sqrt(moy);
}



