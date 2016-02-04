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
  
  Model model;

  model.m=_INDEX_MAX; // num of conservative variables
  model.NumFlux=Gyro_Lagrangian_NumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=Gyro_Lagrangian_BoundaryFlux;
  model.InitData=GyroInitData;
  model.ImposedData=GyroImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;
    
  int deg[]={1, 1, 1};
  int raf[]={1, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);


  InitSimulation(&simu, &mesh, deg, raf, &model);

  simu.vmax = _VMAX; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.2;
  schnaps_real dt = 0;
  schnaps_real tmax = 0.;
  RK2(&simu,tmax);
 
  // save the results and the error
  PlotFields(0,(1==0),&simu,"sol","dgvisu.msh");
  //Plotfield(0,(1==1),&f,"error","dgerror.msh");

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
