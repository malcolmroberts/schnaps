#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "../test.h"
#include "collision.h"
#include "quantities_collision.h"
#include "solverpoisson.h"


void Test_TransportVP_ImposedData(double x[3],double t,double w[]);
void Test_TransportVP_InitData(double x[3],double w[]);
double TransportVP_ImposedKinetic_Data(double x[3],double t,double v);


int main(void) {
  
  // unit tests
    
  int resu=Test_TransportVP();
	 
  if (resu) printf("poisson test OK !\n");
  else printf("poisson test failed !\n");

  return !resu;
} 


int Test_TransportVP(void) {

  bool test=true;

  Field f;

  int vec=1;
  
  f.model.m=_MV+6; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  f.model.vmax = _VMAX; // maximal wave speed
  f.model.NumFlux=Collision_Lagrangian_NumFlux;
  f.model.BoundaryFlux=Collision_Lagrangian_BoundaryFlux;
  f.model.InitData=Test_TransportVP_InitData;
  f.model.ImposedData=Test_TransportVP_ImposedData;
  //f.model.Source = NULL;
  f.model.Source = CollisionSource;
  f.varindex=GenericVarindex;
    
    
  f.interp.interp_param[0]=f.model.m;  // _M
  f.interp.interp_param[1]=3;  // x direction degree
  f.interp.interp_param[2]=0;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=32;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement
 // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  // try to detect a 2d mesh
  bool is1d=Detect1DMacroMesh(&(f.macromesh));
  assert(is1d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  InitField(&f);
  f.macromesh.is1d=true;
  f.is1d=true;

  // prudence...
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  printf("cfl param =%f\n",f.hmin);

  RK2_Poisson(&f,0.03,0.05,0,1,0.0,1.0);

   // save the results and the error
  int iel=_NB_ELEM_V/2;
  int iloc=_DEG_V/2;
  printf("Trace vi=%f\n",-_VMAX+iel*_DV+_DV*glop(_DEG_V,iloc));
  PlotField(iloc+iel*_DEG_V,(1==0),&f,"dgvisu.msh");
  PlotField(iloc+iel*_DEG_V,(1==1),&f,"dgerror.msh");
  

  double dd_Kinetic=L2_Kinetic_error(&f);
  
  printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  test= test && (dd_Kinetic<1e-2);

  return test;

};

void Test_TransportVP_ImposedData(double x[3],double t,double w[]){

  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vi = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j));

    w[i]=TransportVP_ImposedKinetic_Data(x,t,vi);

  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI]=x[0];
  w[_INDEX_EX]=1;
  w[_INDEX_RHO]=0.; //rho init
  w[_INDEX_VELOCITY]=0; // u init
  w[_INDEX_PRESSURE]=0; // p init
  w[_INDEX_TEMP]=0; // e ou T init

};

void Test_TransportVP_InitData(double x[3],double w[]){

  double t=0;
  Test_TransportVP_ImposedData(x,t,w);

};

double TransportVP_ImposedKinetic_Data(double x[3],double t,double v){
  double f;
  double pi=4*atan(1.);
  double xnew=0, vnew=0;
 
  f=exp(-(v-t)*(v-t))*exp(-36*((x[0]-v*t+0.5*t*t)-0.5)*((x[0]-v*t+0.5*t*t)-0.5));
  return f;
};
