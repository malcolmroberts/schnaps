#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "../test/test.h"
#include "quantities_vp.h"
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
  
    
  int deg[]={1, 1, 1};
  int raf[]={1, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  // nbelemv = 10
  // deg_v = 4
  InitKineticData(&schnaps_kinetic_data,10,4);
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  //model.NumFlux=Gyro_Upwind_NumFlux;
  model.NumFlux=GyroZeroNumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=GyroBoundaryFlux;
  model.InitData=GyroInitData;
  model.ImposedData=GyroImposedData;
  //model.Source = NULL;
  model.Source = GyroSource;


  InitSimulation(&simu, &mesh, deg, raf, &model);

  simu.pre_dtfields = UpdateGyroPoisson;
   simu.vmax = kd->vmax; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.2;
  schnaps_real dt = 0;
  schnaps_real tmax = 0;
  RK4(&simu,tmax);
  compute_charge_density(&simu);
  // save the results and the error
  PlotFields(1,(1==0),&simu,"sol","dgvisu.msh");
  PlotFields(kd->index_rho,(1==0),&simu,"sol","dgvisu.msh");
  PlotFields(1,(1==1),&simu,"error","dgerror.msh");

  double dd=L2error(&simu);
  //double dd_l2_vel =GyroL2VelError(&f)
  //double dd_Kinetic=L2_Kinetic_error(&f);
  
  //printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd < 0.005);


  return test; 

};
