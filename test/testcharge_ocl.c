#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "quantities_vp.h"
#include "gyro.h"


int TestChargeOCL(void);
void ChargeOCLInitData(schnaps_real x[3],schnaps_real w[]);
void ChargeOCLImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[]);
void ChargeOCLBoundaryFlux(schnaps_real x[3],schnaps_real t,
		      schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux);
void ChargeOCLZeroNumFlux(schnaps_real wL[],schnaps_real wR[],
			  schnaps_real* vnorm,schnaps_real* flux);
int main(void) {
  
  // unit tests
    
  int resu=TestChargeOCL();
	 
  if (resu) printf("chargeOCL test OK !\n");
  else printf("chargeOCL test failed !\n");

  return !resu;
} 


int TestChargeOCL(void) { 

  bool test=true;

  // read the gmsh file
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../geo/cylindre.msh");

  mesh.period[2]=1;
  BuildConnectivity(&mesh);
    
  int deg[]={2, 2, 1};
  int raf[]={2, 2, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  int nbelemv = 2;
  int deg_v = 2;
  InitKineticData(&schnaps_kinetic_data,nbelemv,deg_v);
  kd->solve_quasineutrality = true;
  kd->substract_mean_charge = false;
  kd->qn_damping = 0;
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  model.NumFlux=ChargeOCLZeroNumFlux;
  //model.NumFlux=GyroZeroNumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=GyroBoundaryFlux;
  model.InitData=ChargeOCLInitData;
  model.ImposedData=ChargeOCLImposedData;
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
  //PlotFields(kd->index_rho,(1==0),&simu,"init_rho","init_rho.msh");
  PlotFields(kd->index_phi,(1==0),&simu,"sol_phi","init_phi.msh");
  //simu.pre_dtfields = UpdateGyroPoisson;
   simu.vmax = kd->vmax; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.2;
  schnaps_real dt = 0.0125;
  schnaps_real tmax = 0.025;
  //RK2(&simu, tmax);
  RK4_CL(&simu,tmax, dt, 0, 0, 0);
  //Computation_charge_density(&simu);
  // save the results and the error
  //PlotFields(1,(1==0),&simu,"sol","dgvisu.msh");
  
  CopyfieldtoCPU(&simu);
  PlotFields(kd->index_phi,(1==0),&simu,"sol_phi","sol_potential.msh");
  PlotFields(kd->index_phi,(1==1),&simu,"err_phi","err_potential.msh"); 

  PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","sol_rho.msh");
  PlotFields(kd->index_rho,(1==1),&simu,"err_rho","err_rho.msh");
  //PlotFields(1,(1==1),&simu,"error","dgerror.msh");

  double dd=L2error(&simu);
  //double dd_l2_vel =ChargeOCLL2VelError(&f)
  //double dd_Kinetic=L2_Kinetic_error(&f);
  
  //printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd < 0.3);


  return test; 

};

void ChargeOCLInitData(schnaps_real x[3],schnaps_real w[]){
  schnaps_real t=0;
  ChargeOCLImposedData(x,t,w);
}

void ChargeOCLImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real pi=4*atan(1.);
  for(int i = 0; i <kd->index_max_kin + 1; i++){
    w[i] =1./12.;
  }
  // exact value of the potential
  // and electric field
  w[kd->index_phi]=-(x[0] * x[0] + x[1] * x[1])/4;
  w[kd->index_rho] =0;
  w[kd->index_ex]=0;//1;
  w[kd->index_ey]=0;//1;
  w[kd->index_ez]=0;
  w[kd->index_u] = 0; 
  w[kd->index_P] = 0; 
  w[kd->index_T] = 0; 

}

void ChargeOCLBoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux)
{
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  ChargeOCLImposedData(x,t,wR);
  GyroUpwindNumFlux(wL,wR,vnorm,flux);
}

void ChargeOCLZeroNumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  
  KineticData *kd = &schnaps_kinetic_data;
  for(int i=0;i<kd->index_max_kin+1;i++){
    flux[i]=0;
  }
  flux[kd->index_phi] = 0;
  flux[kd->index_rho] =0;
  flux[kd->index_ex] = 0;
  flux[kd->index_ey] = 0;
  flux[kd->index_ez] = 0;
  flux[kd->index_u] = 0; 
  flux[kd->index_P] = 0; 
  flux[kd->index_T] = 0;
}
