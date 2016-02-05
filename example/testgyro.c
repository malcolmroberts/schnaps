#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "../test/test.h"
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

  // read the gmsh file
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../geo/cylindre.msh");

  mesh.period[2]=2;
  BuildConnectivity(&mesh);

  int vec=1;
  
    
  int deg[]={3, 3, 3};
  int raf[]={2, 2, 2};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  printf("_MV=%d\n",_MV);
  printf("_INDEX_MAX=%d\n",_INDEX_MAX);
  printf("_INDEX_MAX_KIN=%d\n",_INDEX_MAX_KIN);

  model.m=_INDEX_MAX; // num of conservative variables
  model.NumFlux=Gyro_Upwind_NumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=Gyro_Lagrangian_BoundaryFlux;
  model.InitData=GyroInitData;
  model.ImposedData=GyroImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;


  InitSimulation(&simu, &mesh, deg, raf, &model);

  simu.vmax = _VMAX; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.2;
  schnaps_real dt = 0;
  schnaps_real tmax = 0.1;
  RK4(&simu,tmax);
 
  // save the results and the error
  PlotFields(1,(1==0),&simu,"sol","dgvisu.msh");
  PlotFields(1,(1==1),&simu,"error","dgerror.msh");

  double dd=L2error(&simu);
  //double dd_l2_vel =GyroL2VelError(&f)
  //double dd_Kinetic=L2_Kinetic_error(&f);
  
  //printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  //Velocity_distribution_plot(f.wn);
  test= test && (dd<3e-4);


  return test; 

};
