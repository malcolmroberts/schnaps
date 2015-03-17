#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"



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
  // do not change the potential !
  flux[_MV]=0;
  
};

//! \brief compute square of velocity L2 error
//! \param[in] w : values of f at glops
double L2VelError(double* x,double t,double *w){


  double wex[_MV+1];
  double err2=0;
  CollisionImposedData(x, t,wex);
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      double omega=wglop(_DEG_V,iloc);
      double vi=iel*_DV+_DV*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      err2+=omega*_DV*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
    }
  }
  return err2;
};



void Collision_Lagrangian_BoundaryFlux(double x[3],double t,double wL[],double* vnorm,
				       double* flux){
  double wR[_MV+1];
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
  // exact value of the potential
  w[_MV]=0;
};


double Collision_ImposedKinetic_Data(double x[3],double t,double v){
  double f;
  f=cos(x[0]-v*t);
  return f;
};

double L2_Kinetic_error(field* f){

  double error=0;
  double error_space=0;
  double moy=0; // mean value
  double moy_space=0;


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
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      error+=L2VelError(xphy,f->tnow,w)*wpg*det;
    }
  }
  //moy=moy+weight*moy_space;

  return sqrt(error);
  //return sqrt(error)/sqrt(moy);
}

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void CollisionSource(double* force,double* w, double* source){

  double E=force[0]; // electric field
  double Md[_MV];
  for(int iv=0;iv<_MV;iv++){
    Md[iv]=0;
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
	source[ipg]-=E*omega*w[kpg]*dlag(_DEG_V,iloc,kloc);
      }
    }
  }

  // upwinding
  if (E>0){
    source[_MV-1]+=E*w[_MV-1];
  }
  else {
    source[0]+=-E*w[0];
  }

  for(int iv=0;iv<_MV;iv++){
    source[iv]/=Md[iv];
  }
  

};


void SolvePoisson(field *f){

  // for the moment, works only for the 1d case
  assert(f->macromesh.is1d);

  // assembly of the rigidity matrix

  Skyline sky;

  // number of equation of the Poisson solver
  // = number of nodes in the mesh
  int degx=f->interp.interp_param[1];
  int nelx=f->interp.interp_param[4];
  double xmin=0;
  double xmax=1;  // TO DO: compute the maximal x coordinate
  int neq=degx*nelx+1;
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m=f->model.m;

  InitSkyline(&sky,neq);

  // compute the profile of the matrix
  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	int ino=iloc + ie * degx;
	int jno=jloc + ie * degx;
	SwitchOn(&sky,ino,jno);
      }
    }
  }

  AllocateSkyline(&sky);

  // local matrix (assuming identical finite elements)
  double aloc[degx+1][degx+1];
  for(int iloc=0;iloc<degx+1;iloc++){
    for(int jloc=0;jloc<degx+1;jloc++){
      aloc[iloc][jloc]=0;
    }
  }
  for(int ipg=0;ipg<degx+1;ipg++){
    double omega=wglop(degx,ipg);
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	double dxi=dlag(degx,iloc,ipg);
	double dxj=dlag(degx,jloc,ipg);
	aloc[iloc][jloc]+=dxi*dxj*omega*nelx;
      }
    }
  }

  // assembly of the matrix
  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	int ino=iloc + ie * degx;
	int jno=jloc + ie * degx;
	double val = aloc[iloc][jloc];
	SetSkyline(&sky,ino,jno,val);
      }
    }
  }

  // source assembly 
  double source[neq];
  for(int i=0;i<neq;i++){
    source[i]=0;
  }
  double cst=1e0; // constant source TO DO: replace by the charge
  for(int ie=0;ie<nelx;ie++){
    for(int ipg=0;ipg<degx+1;ipg++){
      double omega=wglop(degx,ipg);
      for(int iloc=0;iloc<degx+1;iloc++){
	int ino=iloc + ie * degx;
	source[ino]+= cst*omega/nelx;
      }
    }
  }



  // dirichlet boundary condition at the first and last location
  SetSkyline(&sky,0,0,1e20);
  SetSkyline(&sky,neq-1,neq-1,1e20);

  //DisplaySkyline(&sky);

  FactoLU(&sky);

  double sol[neq];
  SolveSkyline(&sky,source,sol);

  for(int i=0;i<neq;i++){
    printf("sol %d = %f\n",i,sol[i]);
  }

  // now put the solution at the right place
  for(int ie=0;ie<nelx;ie++){
     for(int ipg=0;ipg<degx+1;ipg++){
       // position in the continuous vector
       int ino=ipg + ie * degx;
       // position in the DG vector
       int imem=f->varindex(f->interp_param,0,ipg+ie*(degx+1),m-1);
       f->wn[imem]=sol[ino];
     }
  }
	
  FreeSkyline(&sky);
}

