#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
//#include "collision.h"
#include "gyro.h"


int main(void) {
  
  // unit tests
    
  int resu=TestGyro();
	 
  if (resu) printf("gyro test OK !\n");
  else printf("gyro test failed !\n");

  return !resu;
} 


int TestGyro(void) {

  bool test=true;

  Field f;

  int vec=1;
  
  f.model.m=_MV; // num of conservative variables
  f.model.vmax = _VMAX; // maximal wave speed 
  f.model.NumFlux=Gyro_Lagrangian_NumFlux;
  f.model.BoundaryFlux=Gyro_Lagrangian_BoundaryFlux;
  f.model.InitData=GyroInitData;
  f.model.ImposedData=GyroImposedData;
  f.model.Source = NULL;
  f.varindex=GenericVarindex;
    
    
  f.interp.interp_param[0]=_MV;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=2;  // z direction degree
  f.interp.interp_param[4]=4;  // x direction refinement
  f.interp.interp_param[5]=4;  // y direction refinement
  f.interp.interp_param[6]=4;  // z direction refinement
  // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  // try to detect a 2d mesh
  //bool is1d=Detect1DMacroMesh(&(f.macromesh));
  //assert(is1d);

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

  // time derivative
  //dtField(&f);
  //DisplayField(&f);
  //assert(1==2);
  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  RK2(&f,0.1);
 
  // save the results and the error
  PlotField(0,(1==0),&f,"dgvisu.msh");
  PlotField(0,(1==1),&f,"dgerror.msh");

  double dd=L2error(&f);
  double dd_Kinetic=GyroL2_Kinetic_error(&f);
  
  printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  test= test && (dd<3e-4);


  return test;

};


