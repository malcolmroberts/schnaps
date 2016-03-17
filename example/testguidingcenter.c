#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "quantities_vp.h"
#include "gyro.h"
#include "solverpoisson.h"
#include "diagnostics_vp.h"

int TestGuidingCenter(void);

void SolveQuasineutreEq(void *field);
void GuidingCInitData(schnaps_real x[3],schnaps_real w[]);
void GuidingCImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[]);
void GuidingCBoundaryFlux(schnaps_real x[3],schnaps_real t,
		      schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux);
//void SolveQuasineutreEq(void *field, schnaps_real *w, schnaps_real dt); 
int main(void) {
  
  // unit tests
    
  int resu=TestGuidingCenter();
	 
  if (resu) printf("guiding center test OK !\n");
  else printf("guiding center test failed !\n");

  return !resu;
} 


int TestGuidingCenter(void) { 

  bool test=true;

  // read the gmsh file
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../geo/cylindre.msh"); 
  //ReadMacroMesh(&mesh,"testcylindre.msh");
  /* Detect2DMacroMesh(&mesh); */
  /* bool is2d=mesh.is2d; */
  /* assert(is2d); */
  mesh.period[2]=2;
  BuildConnectivity(&mesh);
  
    
  int deg[]={3, 3, 1};
  int raf[]={5, 5, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  int nbelemv = 1;
  int deg_v = 0;
  
  InitKineticData(&schnaps_kinetic_data,nbelemv,deg_v);
  /* kd->vmax = 1; */
  /* kd->dv = 1; */
  kd->solve_quasineutrality = true;
  kd->substract_mean_charge = false;
  kd->qn_damping = 0;
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  model.NumFlux=GyroUpwindNumFlux;
  //model.NumFlux=GyroZeroNumFlux;
  model.BoundaryFlux=GuidingCBoundaryFlux;
  model.InitData=GuidingCInitData;
  model.ImposedData=GuidingCImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;


  InitSimulation(&simu, &mesh, deg, raf, &model);
  //check the initial distribution funtion
  PlotFields(kd->index_rho,(1==0),&simu,"init_rho","testres2/initrho.msh");
  PlotFields(kd->index_max_kin,(1==0),&simu,"init_f","testres2/initdistrib.msh");
   simu.vmax = kd->vmax; // maximal wave speed
   simu.pre_dtfields = SolveQuasineutreEq;
   /* simu.post_dtfields=NULL; */
  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.7;
  //schnaps_real dt = 0;
  schnaps_real tmax = 10;
  RK2(&simu,tmax);
  PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","testres3/rho.msh");
  PlotFields(kd->index_max_kin,(1==0),&simu,"sol_f","testres3/distrib.msh");
  PlotFields(kd->index_ex,(1==0),&simu,"sol_ex","testres3/ex.msh");
  PlotFields(kd->index_ey,(1==0),&simu,"sol_ey","testres3/ey.msh");
  PlotFields(kd->index_ez,(1==0),&simu,"sol_ez","testres3/ez.msh");
  PlotFields(kd->index_phi,(1==0),&simu,"sol_phi","testres3/potential.msh");
  //Plot_Energies(&simu, simu.dt);
  /* PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","testres2/rho.msh"); */
  /* PlotFields(kd->index_max_kin,(1==0),&simu,"sol_f","testres2/distrib.msh"); */
  /* PlotFields(kd->index_ex,(1==0),&simu,"sol_ex","testres2/ex.msh"); */
  /* PlotFields(kd->index_ey,(1==0),&simu,"sol_ey","testres2/ey.msh"); */
  /* PlotFields(kd->index_ez,(1==0),&simu,"sol_ez","testres2/ez.msh"); */
  /* PlotFields(kd->index_phi,(1==0),&simu,"sol_phi","testres2/potential.msh"); */
  /* PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","rho.msh"); */
  /* PlotFields(kd->index_max_kin,(1==0),&simu,"sol_f","distrib.msh"); */
  /* PlotFields(kd->index_ex,(1==0),&simu,"sol_ex","ex.msh"); */
  /* PlotFields(kd->index_ey,(1==0),&simu,"sol_ey","ey.msh"); */
  /* PlotFields(kd->index_ez,(1==0),&simu,"sol_ez","ez.msh"); */
  /* PlotFields(kd->index_phi,(1==0),&simu,"sol_phi","potential.msh"); */
  //Plot_Energies(&simu, simu.dt);
  double dd=L2error(&simu);
  //double dd_l2_vel =GyroL2VelError(&f)
  double dd_Kinetic=L2_Kinetic_error(&simu);
  
  printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd_Kinetic < 0.005);


  return test; 

};


void GuidingCInitData(schnaps_real x[3],schnaps_real w[]){
  schnaps_real t=0;
  GuidingCImposedData(x,t,w);
}
void GuidingCImposedData1(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData *kd = &schnaps_kinetic_data;
  //For the same with drif kinetic
  schnaps_real m=5;
  schnaps_real eps = 1e-6;
  schnaps_real pi= 4.0*atan(1.0);
  schnaps_real r = sqrt(x[0] * x[0] + x[1] * x[1]);
  schnaps_real phi = atan(x[1]/x[0]);
 
  /* if (x[0]*x[1] < 0){ */
  /*   phi += pi; */
  /*  } */
  for(int i = 0; i <kd->index_max_kin + 1; i++){
       w[i] = 1+eps*exp(-(r-7.3)*(r-7.3)/8)*cos(m*phi);
  }
 // exact value of the potential
  // and electric field
  w[kd->index_phi]=0; //(x[0] * x[0] + x[1] * x[1])/4;
  w[kd->index_rho] = 0; //-1 + kd->qn_damping * w[kd->index_phi];
  w[kd->index_rho] = 1+eps*exp(-(r-7.3)*(r-7.3)/8)*cos(m*phi);//}
  w[kd->index_ex]=0; //-x[0]/2;
  w[kd->index_ey]=0; //-x[1]/2;
  w[kd->index_ez]=0;
  w[kd->index_u] = 0; // u init
  w[kd->index_P] = 0; // p init
  w[kd->index_T] = 0; // e ou T init
}


void GuidingCImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData *kd = &schnaps_kinetic_data;
  
  //anneaux
  schnaps_real m=6;
  schnaps_real eps = 0.001;
  schnaps_real pi= 4.0*atan(1.0);
  schnaps_real r = sqrt(x[0] * x[0] + x[1] * x[1]);
  schnaps_real phi = atan(x[1]/x[0]);
  schnaps_real rmin = 5;
  schnaps_real rmax = 8;
  /* if (x[0]*x[1] < 0){ */
  /*   phi += pi; */
  /* } */
 for(int i = 0; i <kd->index_max_kin + 1; i++){
      w[i] =exp(-(r-6.5)*(r-6.5)*4)*(1+eps*cos(m*phi));
  }

  w[kd->index_rho] = exp(-(r-6.5)*(r-6.5)*4)*(1+eps*cos(m*phi));
  /////////////////////////////////////////////
 /*  schnaps_real rmin = 4; */
 /*  schnaps_real rmax = 6; */
 /*  if (x[0]*x[1] < 0){ */
 /*    phi += pi; */
 /*  } */
 /* for(int i = 0; i <kd->index_max_kin + 1; i++){ */
 /*      w[i] =exp(-(r-5)*(r-5)*4)*(1+eps*cos(m*phi)); */
 /*  } */
 /* //w[kd->index_rho] = exp(-(r-5)*(r-5)*4)*(1+eps*cos(m*phi)); */
   
 /*  /\* for(int i = 0; i <kd->index_max_kin + 1; i++){ *\/ */
 /*  /\*   w[i] =exp(-(r-6.5)*(r-6.5)*4)*(1+eps*cos(m*phi)); *\/ */
 /*  /\* } *\/ */
 /*  w[kd->index_rho] = exp(-(r-6.5)*(r-6.5)*4)*(1+eps*cos(m*phi));  */
  w[kd->index_phi]=0;
  w[kd->index_ex]=0; 
  w[kd->index_ey]=0;
  w[kd->index_ez]=0;
  w[kd->index_u] = 0; // u init
  w[kd->index_P] = 0; // p init
  w[kd->index_T] = 0; // e ou T init
}


void GuidingCBoundaryFlux(schnaps_real x[3],schnaps_real t,
		      schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux)
{
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  GuidingCImposedData(x,t,wR);
  GyroUpwindNumFlux(wL,wR,vnorm,flux);
}


void SolveQuasineutreEq(void *si) {
  Simulation *simu = si;
  KineticData * kd=&schnaps_kinetic_data;
  int type_bc = 1;

  Computation_charge_density(simu);
  static ContinuousSolver ps;
  static bool is_init = false;

  if (!is_init){
    is_init = true;
    int nb_var=1;
    int * listvar= malloc(nb_var * sizeof(int));
    listvar[0]=kd->index_phi;
    InitContinuousSolver(&ps,simu,type_bc,nb_var,listvar);

    ps.matrix_assembly=ContinuousOperator_Poisson2D;
    ps.rhs_assembly=RHSPoisson_Continuous;
    ps.bc_assembly= ExactDirichletContinuousMatrix;
    ps.postcomputation_assembly=Computation_ElectricField_Poisson;
    //ps.postcomputation_assembly=NULL;

#undef PARALUTION
#ifdef PARALUTION
    ps.lsol.solver_type = PAR_LU;
    ps.lsol.pc_type=NONE;
#else
    ps.lsol.solver_type = LU;
    ps.lsol.pc_type=NONE;
#endif
  }
 
  SolveContinuous2D(&ps);
  //freeContinuousSolver(&ps);
}
  
 
