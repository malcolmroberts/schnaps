//#include "collision.h"
#include "gyro.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"


/* void Gyro_Lagrangian_NumFlux(double wL[],double wR[],double* vnorm,double* flux){ */
  
/*   for(int i=0;i<_MV;i++){ */
/*     int j=i%_DEG_V; // local connectivity put in function */
/*     int nel=i/_DEG_V; // element num (TODO : function) */

/*     double vn = (nel*_DV + */
/* 		 _DV* glop(_DEG_V,j))*vnorm[0]; */
    
/*     double vnp = vn>0 ? vn : 0; */
/*     double vnm = vn-vnp; */

/*     flux[i] = vnp * wL[i] + vnm * wR[i]; */
/*   } */
  
/* }; */

//flux num gyro
void Gyro_Lagrangian_NumFlux(double wL[],double wR[],double* vnorm,double* flux){
  double eps =0; //if not equal 0 => decentered flux
  /* double E_x =0; //firstly consider the electric field is const */
  /* double E_y =1; */
  double E_x =wL[_INDEX_EX];
  double E_y =wL[_INDEX_EY];
  /* printf("Ey = %f \n",E_y); */
  /* assert(1==2); */
  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)  
    double wm = (wL[i]+wR[i])/2;
    double flux1 = E_y*wm;
    double flux2 = -E_x*wm;
    double v =-_VMAX+ nel*_DV +
      _DV* glop(_DEG_V,j); // gauss_lob_point[j]
    double flux3 =v*wm;
  
    flux[i] = vnorm[0]*flux1+vnorm[1]*flux2+vnorm[2]*flux3-eps*(wR[i]-wL[i])/2;
  }
  flux[_INDEX_PHI] =0;
  flux[_INDEX_EX]=0;
  flux[_INDEX_EY]=0;
  
};

//! \brief compute square of velocity L2 error
//! \param[in] w : values of f at glops
double GyroL2VelError(double* x,double t,double *w){


  double wex[_MV];
  double err2=0;
  GyroImposedData(x, t,wex);
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



void Gyro_Lagrangian_BoundaryFlux(double x[3],double t,double wL[],double* vnorm,
				       double* flux){
  double wR[_MV+3];
  GyroImposedData(x,t,wR);
  Gyro_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};


void GyroInitData(double x[3],double w[]){

  double t=0;
  GyroImposedData(x,t,w);

};



void GyroImposedData(double x[3],double t,double w[]){

  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vi = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j));
    //w[i]=cos(x[0]-vi*t);
    double pi=4*atan(1.);
    double xi = x[0]-t;
    //pour transport 1d
    // w[i] = cos(2*pi*xi);//*exp(-0.5*pow(xi-0.5,2));
    //w[i] = cos(2*pi*xi)*exp(-1/(1-pow(2*xi-1,2)));
    //w[i] = exp(-4*pow(xi-0.5,2));
    double yi = x[1]+t;
    //pour transport 2D
    //w[i] = exp(-4*pow(yi-0.5,2))*exp(-4*pow(xi-0.5,2));
    //w[i] = cos(2*pi*yi);//*exp(-1/(1-pow(2*xi-1,2)));
    //w[i] = cos(2*pi*xi)*cos(2*pi*yi);//*exp(-0.5*pow(xi-0.5,2));
        //pour transport 3D
     double zi=x[2]-vi*t;
     w[i] = exp(-4*pow(yi-0.5,2))*exp(-4*pow(xi-0.5,2))*exp(-4*pow(zi-0.5,2));
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI]=0;
  w[_INDEX_EX]=1;//1;
  w[_INDEX_EY]=1;//1;

};


double Gyro_ImposedKinetic_Data(double x[3],double t,double v){
  double f;
  f=cos(x[0]-v*t);
  return f;
};

double GyroL2_Kinetic_error(Field* f){

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
      error+=GyroL2VelError(xphy,f->tnow,w)*wpg*det;
    }
  }
  //moy=moy+weight*moy_space;

  return sqrt(error);
  //return sqrt(error)/sqrt(moy);
}

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void GyroSource(double* force,double* w, double* source){

  double E=force[0]; // electric field
  double Md[_MV];
  for(int iv=0;iv<_MV;iv++){
    Md[iv]=0;
    source[iv]=0;
  }
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


