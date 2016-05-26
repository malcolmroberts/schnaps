#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "quantities_vp.h"
#include "gyro.h"


int TestGyro(void);

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

  mesh.period[2]=1;
  BuildConnectivity(&mesh);

  
    
  int deg[]={4, 4, 2};
  int raf[]={3, 3, 2};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  int nbelemv = 4;
  int deg_v = 2;
  InitKineticData(&schnaps_kinetic_data,nbelemv,deg_v);
  kd->solve_quasineutrality = true;
  kd->substract_mean_charge = false;
  kd->qn_damping = 0;
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  //model.NumFlux=GyroUpwindNumFlux;
  model.NumFlux=GyroZeroNumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=GyroBoundaryFlux;
  model.InitData=GyroInitData;
  model.ImposedData=GyroImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;

  
  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _NBELEMV=%d ",nbelemv);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _DEGV=%d ",deg_v);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _VMAX=%f ", kd->vmax);
  strcat(cl_buildoptions, buf);


  schnaps_ocl_getcharge = true;
  InitSimulation(&simu, &mesh, deg, raf, &model);
  PlotFields(kd->index_phi,(1==0),&simu,"init_phi","init_potential.msh");
  //simu.pre_dtfields = UpdateGyroPoisson;
   simu.vmax = kd->vmax; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.2;
  schnaps_real dt = 0;
  schnaps_real tmax = 0.0001;
  //RK2(&simu,tmax);
  RK4_CL(&simu,tmax, dt, 0, 0, 0);
  //Computation_charge_density(&simu);
  // save the results and the error
  //PlotFields(1,(1==0),&simu,"sol","dgvisu.msh");
  
  CopyfieldtoCPU(&simu);
  PlotFields(kd->index_phi,(1==0),&simu,"sol","dgvisu.msh");
  PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","sol_rho.msh");
  PlotFields(kd->index_ex,(1==0),&simu,"sol_ex","sol_ex.msh");
  PlotFields(kd->index_ey,(1==0),&simu,"sol_ey","sol_ey.msh");
  PlotFields(kd->index_ez,(1==0),&simu,"sol_ez","sol_ez.msh");
  //PlotFields(1,(1==1),&simu,"error","dgerror.msh");

  double dd=L2error(&simu);
  //double dd_l2_vel =GyroL2VelError(&f)
  //double dd_Kinetic=L2_Kinetic_error(&simu);
  
  
  //printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd < 0.005);


  return test; 

};
