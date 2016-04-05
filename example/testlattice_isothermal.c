#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "global.h"
#include "lattice.h"
#include "implicit.h"
//#include "field.h"
//
typedef struct LbmSimuParams{
  int degx;
  int degy;
  int rafx;
  int rafy;
  schnaps_real cfl;
  schnaps_real tmax;
  schnaps_real cref;
  schnaps_real tau;
  schnaps_real diag_2d_period;
} LbmSimuParams;
typedef struct DoubleShearKHParams{
  schnaps_real kappa;
  schnaps_real delta;
  schnaps_real uref;
} DoubleShearKHParams;
typedef struct Linear2DWaveParams{
  int nkx;
  int nky;
  schnaps_real offset;
} Linear2DWaveParams;
//
int TestLattice_isothermal_DoubleShearKH(void);
void DoubleShearKH_InitData(schnaps_real x[3], schnaps_real w[]);
void DoubleshearKH_Plot_Fields(void *s,schnaps_real *w);
//
int TestLattice_isothermal_Linear2DWave(void);
int TestLattice_isothermal_Linear2DWaveImplicit(void);
void Linear2DWave_InitData(schnaps_real x[3],schnaps_real w[]);
void Linear2DWave_ImposedData(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]);
void Linear2DWave_ImposedData_OneNode(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]);
void Linear2DWave_Periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
               schnaps_real *flux);
void Linear2DWave_Periodic_BoundaryFlux_OneNode(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
               schnaps_real *flux);
void Linear2DWave_Plot_Fields(void *s,schnaps_real *w);
void Linear2DWave_CollectDiags(void *simu,schnaps_real *diag_vals);
//
void Relaxation(void* s);
void Moments(void * s);
void Vorticity_2D_computation(Simulation *simu, int ifield_ux, int ifield_uy);
void Compute_and_dump_Vorticity_2d(Simulation *simu,char *filename, int create_file, schnaps_real t, int istep);
//
void PlotExtFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
void PlotExtFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
void PlotFieldsBinSparseMultitime(int typplot, int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep);
void PlotExtScalarFieldBinMultitime(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename,int create_file, schnaps_real t,int istep);
void PlotVecFieldsBinSparseMultitime(int typplot[3], int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep);
//
void Store_Lattice_diags(Simulation *simu);
void Dump_Lattice_Diagnostics(Simulation *simu,char simtag[3]);
//
void Lattice_Dummy_m1_InitData(schnaps_real x[3],schnaps_real w[]);
void LatticeThetaTimeScheme_MultiSim(Simulation *simu, Simulation **simu_advec,schnaps_real tmax, schnaps_real dt);
void LatticeThetaTimeScheme(Simulation *simu,schnaps_real tmax, schnaps_real dt);

// global parameters with default values
LbmSimuParams SimParams={
  .degx=4,.degy=4,
  .rafx=8,.rafy=8,
  .cfl=1.0,.tmax=1.0,
  .cref=1.0,.tau=0.0001,
  .diag_2d_period=1.0};
DoubleShearKHParams DKHParams={.kappa=80.0,.delta=0.05,.uref=0.05};
Linear2DWaveParams  LW2DParams={.nkx=1,.nky=0,.offset=0.0};
//
// global dat for taggig diags (this should be in a structure somewhere simu ?)
char simutag[3]="TAG";
//
int main(void) {
#define LBMDOUBLESHEARKH 0
#define LBMWAVE2D 1
//
#define LBMTESTCASE 1
#if (LBMTESTCASE==LBMWAVE2D)
  printf(" LBM - 2D - Isothermal - Linear2DWave\n");
  //
  SimParams.degx=4;
  SimParams.degy=4;
  SimParams.rafx=8;
  SimParams.rafy=8;
  
  SimParams.cfl=1.0;
  SimParams.tmax=0.05;
  SimParams.tau=0.000000001;
  SimParams.diag_2d_period=0.1;
  //
  LW2DParams.nkx=3;
  LW2DParams.nky=2;
  LW2DParams.offset=0.0;
  // 
  //int resu= TestLattice_isothermal_Linear2DWave();
  int resu= TestLattice_isothermal_Linear2DWaveImplicit();
#endif
#if (LBMTESTCASE==LBMDOUBLESHEARKH)
  printf(" LBM - 2D - Isothermal - Double sheared flow  KH\n"); 
  int resu=TestLattice_isothermal_DoubleShearKH();
#endif
  if (resu) printf("lattice test OK !\n");
  else printf("lattice test failed !\n");
  return !resu;
} 
//
// Basic routines common to all tests
void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  
  LatticeData * ld=&schnaps_lattice_data;

};


void Relaxation(void* s){
  Simulation * simu =s;
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real w_eq[simu->wsize];
  Compute_distribution_eq(simu,w_eq);
  Compute_relaxation(simu,w_eq);
}

void Moments(void* s){
  Simulation * simu =s;
  
  Compute_moments(simu);
}
/*************************************************************************************************/
/* KH Double shear test routines*/
/*************************************************************************************************/
int TestLattice_isothermal_DoubleShearKH(void) {
  
  int test=0;
  //int vec=1;
  //schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);

  
  Model model;
  LatticeData * ld=&schnaps_lattice_data;  
  ld->feq=&feq_isothermal_D2Q9;
  schnaps_real cref= SimParams.cref;
  InitLatticeData(ld,2,9,0,cref);
  ld->diag_2d_period=SimParams.diag_2d_period;
  ld->tau=SimParams.tau;
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  //
  model.NumFlux=Lattice_NumFlux;
  //
  model.InitData = DoubleShearKH_InitData;
  model.ImposedData = NULL;
  model.BoundaryFlux = NULL;
  model.Source = NULL;
  
  int deg[]={SimParams.degx, SimParams.degy, 0};
  int raf[]={SimParams.rafx,SimParams.rafy, 1};
  //
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  //
  InitSimulation(&simu, &mesh, deg, raf, &model);
  Moments(&simu);
  //simu.vmax = 2*ld->c; 
  simu.vmax = ld->c *sqrt(2.0); 
  simu.cfl=SimParams.cfl;
  simu.nb_diags = 1;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk = DoubleshearKH_Plot_Fields;
  schnaps_real tmax = SimParams.tmax;
  //
  DoubleshearKH_Plot_Fields(&simu, NULL);
  //
  RK2(&simu, tmax);
  //
  test= 1;
  return test; 
}
/*KH Double shear Init*/
void DoubleShearKH_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real delta = DKHParams.delta, kappa=DKHParams.kappa,uref=DKHParams.uref;
  //
  schnaps_real rho = 1.0;
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld->c)/3.0;
  schnaps_real p =1.0;
  //
  //printf(" c0 :%f \t temp:%f\n",ld->c, temp);
  if(x[1]<0.5){
    ux=uref * tanh(kappa*(x[1]-0.25));
  }
  else{
    ux=uref * tanh(kappa*(0.75-x[1]));
  }    
  uy = uref * delta * sin(2.0 * my_pi*(x[0]+0.25));
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
};
/********** KH Double shear data dumps
/**************************************************************************************/
void DoubleshearKH_Plot_Fields(void *s,schnaps_real *w){
  Simulation *simu=s;
  LatticeData *ld=&schnaps_lattice_data;
  schnaps_real period=ld->diag_2d_period;
  schnaps_real dt=simu->dt;
  int diagperiod = (int) (period/dt);
  int istep=simu->iter_time_rk;
  schnaps_real t=simu->tnow;
  int tmax=simu->tmax;
  int create_file=0;
  if (istep==0){
    create_file=1;
  }
  else
  {
  create_file = 0;
  }
  if (diagperiod || create_file){
  if ((istep%diagperiod ==0) || (t== tmax)){
    printf("Dumping fields at it=%i\n",istep);
    PlotFieldsBinSparseMultitime(ld->index_rho,false,simu,"rho","lbm_KH2shear_rho.msh",create_file,t,istep);
    PlotFieldsBinSparseMultitime(ld->index_ux,false,simu,"ux","lbm_KH2shear_ux.msh",create_file,t,istep);
    PlotFieldsBinSparseMultitime(ld->index_uy,false,simu,"uy","lbm_KH2shear_uy.msh",create_file,t,istep);
    Compute_and_dump_Vorticity_2d(simu,"lbm_KH2shear_uvort.msh",create_file,t,istep);
  };
  };
};
/**************************************************************************************/
/***************************** 2D Wave equation linear LBM ***************************/ 
int TestLattice_isothermal_Linear2DWave(void) {
  
  int test=0;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);
  //
  Model model;
  LatticeData * ld=&schnaps_lattice_data;  
  ld->feq=&feq_isothermal_linearwave_D2Q9;
  ld->collect_diags=&Linear2DWave_CollectDiags;
  //
  //schnaps_real csound= sqrt(3.0/2.0);
  InitLatticeData(ld,2,9,0,SimParams.cref);
  ld->diag_2d_period=SimParams.diag_2d_period;
  ld->tau=SimParams.tau;
  //  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  //
  model.NumFlux=Lattice_NumFlux;
  //
  model.InitData = Linear2DWave_InitData;
  model.ImposedData = Linear2DWave_ImposedData;
  //model.BoundaryFlux = Linear2DWave_Periodic_BoundaryFlux;
  model.BoundaryFlux =NULL;
  model.Source = NULL;
  //
  int deg[]={SimParams.degx,SimParams.degy, 0};
  int raf[]={SimParams.rafx,SimParams.rafy,1};
  //
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  //
  InitSimulation(&simu, &mesh, deg, raf, &model);
  Moments(&simu);
  simu.vmax = ld->c *sqrt(2.0); 
  simu.cfl=SimParams.cfl;
  simu.nb_diags = 3;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk = Linear2DWave_Plot_Fields;
  schnaps_real tmax = SimParams.tmax;
  //
  //Linear2DWave_Plot_Fields(&simu,NULL);
  sprintf(simutag,"RK2");
  RK2(&simu, tmax);
  //
  schnaps_real end_diags[simu.nb_diags];
  Linear2DWave_CollectDiags(&simu,end_diags);
  printf("RK2 end error:t=%f \tabs : %f \tmean : %f \trel:%f\n",simu.tnow,end_diags[0],end_diags[1],end_diags[2]);
  //
  Dump_Lattice_Diagnostics(&simu,simutag);
  test= 1;
  return test; 
}
int TestLattice_isothermal_Linear2DWaveImplicit(void) {
  
  int test=0;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);
  //
  Model model;
  LatticeData * ld=&schnaps_lattice_data;  
  ld->feq=&feq_isothermal_linearwave_D2Q9;
  ld->collect_diags=&Linear2DWave_CollectDiags;
  //
  //schnaps_real csound= sqrt(3.0/2.0);
  InitLatticeData(ld,2,9,0,SimParams.cref);
  ld->diag_2d_period=SimParams.diag_2d_period;
  ld->tau=SimParams.tau;
  //  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  //
  model.NumFlux=Lattice_NumFlux;
  //
  model.InitData = Linear2DWave_InitData;
  model.ImposedData = Linear2DWave_ImposedData;
  model.BoundaryFlux = Linear2DWave_Periodic_BoundaryFlux;
  model.Source = NULL;
  //
  int deg[]={SimParams.degx,SimParams.degy, 0};
  int raf[]={SimParams.rafx,SimParams.rafy,1};
  //
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  //
  InitSimulation(&simu, &mesh, deg, raf, &model);
  Moments(&simu);
  simu.vmax = ld->c *sqrt(2.0); 
  simu.cfl=SimParams.cfl;
  simu.nb_diags = 3;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk=Linear2DWave_Plot_Fields;
  schnaps_real tmax = SimParams.tmax;
  //
  // basic explicit RK2
  sprintf(simutag,"RK2");
  RK2(&simu, tmax);
  schnaps_real end_diags[simu.nb_diags];
  Linear2DWave_CollectDiags(&simu,end_diags);
  printf("RK2 end error: abs : %f \tmean : %f \trel:%f\n",end_diags[0],end_diags[1],end_diags[2]);
  Dump_Lattice_Diagnostics(&simu,simutag);
  //
  Model model_advec;
  model_advec.m=1;
  model_advec.InitData =Lattice_Dummy_m1_InitData;
  model_advec.NumFlux=Lattice_OneNodeNumFlux;
  model_advec.ImposedData = NULL;
  model_advec.BoundaryFlux = Linear2DWave_Periodic_BoundaryFlux_OneNode;
  model_advec.Source = NULL;
  //
  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  Moments(&simu2);
  simu2.vmax = ld->c *sqrt(2.0); 
  simu2.cfl=SimParams.cfl; // no impact we impose dt
  simu2.nb_diags = 3;
  simu2.pre_dtfields = Relaxation;
  simu2.post_dtfields = Moments;
  simu2.update_after_rk=Linear2DWave_Plot_Fields;
  //simu2.update_after_rk=NULL;
  tmax = SimParams.tmax;
  schnaps_real dt= simu.dt; // we use the timestep of previous RK2 run
  printf("Testing Implicit scheme with dt=%f\n",dt);
  //LatticeThetaTimeScheme_MultiSim(&simu2,simu_advection,tmax,dt);
  sprintf(simutag,"IMP");
  LatticeThetaTimeScheme(&simu2,tmax,dt);
  //*****************************************************************/
  // test implicit
  schnaps_real end_diags2[simu2.nb_diags];
  Linear2DWave_CollectDiags(&simu2,end_diags2);
  //
  printf("Summary\n");
  printf("RK2  :t=%f \t abs : %f \tmean : %f \trel:%f\n",simu.tnow,end_diags[0],end_diags[1],end_diags[2]);
  printf("Imp  :t=%f \t abs : %f \tmean : %f \trel:%f\n",simu2.tnow,end_diags2[0],end_diags2[1],end_diags2[2]);
  //
  Dump_Lattice_Diagnostics(&simu2,simutag);
  //
  test= 1;
  return test; 
}
void Lattice_Dummy_m1_InitData(schnaps_real x[3],schnaps_real w[]){
  w[0]=0.0;
}
void Linear2DWave_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  // wave mode numbers in half integer units
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset=LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky= 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix= kx * x[0];
  schnaps_real phiy= ky * x[1];
  //
  schnaps_real rho = offset+cos(phix) * cos(phiy);
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
}
void Linear2DWave_ImposedData(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]){
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset=LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky= 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix= kx * x[0];
  schnaps_real phiy= ky * x[1];
  //
  schnaps_real k = sqrt(kx * kx + ky * ky);
  schnaps_real c0 = ld->c * sqrt(1/3.0);
  schnaps_real omega= k* c0;
  schnaps_real phit= omega *t;
  //
  //schnaps_real rho = 1.0 + pert * cos(phix);
  schnaps_real rho = offset+cos(phix) * cos(phiy) * cos(phit);
  //schnaps_real ux =uref * pert * sin(phiy);
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
}
void Linear2DWave_ImposedData_OneNode(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]){
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset=LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky= 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix= kx * x[0];
  schnaps_real phiy= ky * x[1];
  //
  schnaps_real k = sqrt(kx * kx + ky * ky);
  schnaps_real c0 = ld->c * sqrt(1/3.0);
  schnaps_real omega= k* c0;
  schnaps_real phit= omega *t;
  //
  //schnaps_real rho = 1.0 + pert * cos(phix);
  schnaps_real rho = offset+cos(phix) * cos(phiy) * cos(phit);
  //schnaps_real ux =uref * pert * sin(phiy);
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  int i_node= ld->current_node_index;
  w[0]= ld->feq(i_node,ld,rho,ux,uy,uz,temp,p);
}
/***************************************************************************************/
void Linear2DWave_Periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
	LatticeData *ld=&schnaps_lattice_data;
  schnaps_real wR[ld->index_max];
  Linear2DWave_ImposedData(x,t,wR);
  Lattice_NumFlux(wL,wR,vnorm,flux);
}
void Linear2DWave_Periodic_BoundaryFlux_OneNode(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux){ 
	LatticeData *ld=&schnaps_lattice_data;
	int i_node=ld->current_node_index;
  schnaps_real wR[1];
  Linear2DWave_ImposedData_OneNode(x,t,wR);
  Lattice_OneNodeNumFlux(wL,wR,vnorm,flux);
}
/**************************************************************************************/
void Linear2DWave_Plot_Fields(void *s,schnaps_real *w){
  Simulation *simu=s;
  LatticeData *ld=&schnaps_lattice_data;
  schnaps_real period=ld->diag_2d_period;
  schnaps_real dt=simu->dt;
  int diagperiod = (int) (period/dt);
  int istep=simu->iter_time_rk;
  schnaps_real t=simu->tnow;
  int tmax=simu->tmax;
  int create_file=0;
  Store_Lattice_diags(simu);
  //printf("Called with istep=%i t=%f cfl=%f\n",istep,t,simu->cfl);
  if (istep==0){
    create_file=1;
  }
  else
  {
  create_file = 0;
  }
  if (diagperiod || create_file){
  if ((istep%diagperiod ==0) || (t== tmax)){
    istep=istep+1;
    printf("Dumping fields at it=%i (period %i)\n",istep,diagperiod);
    int raf=simu->fd[0].raf[0];
    schnaps_real cfl=simu->cfl;
    char filename_rho[sizeof("lbm2DWave_rho_TAG_raf000_cfl0.000.msh")];
    sprintf(filename_rho,"lbm_2DWave_rho_%s_raf%03d_cfl%1.3f.msh",simutag,raf,cfl);
    //char filename_rho_error[sizeof("lbm2DWave_rho_000.msh")];
    //sprintf(filename_rho_error,"lbm_2DWave_rho_error_%03d.msh",raf);
    PlotFieldsBinSparseMultitime(ld->index_rho,false,simu,"rho",filename_rho,create_file,t,istep);
    PlotFieldsBinSparseMultitime(ld->index_rho,true,simu,"rho_error",filename_rho,0,t,istep);
  };
  }
};
/**************************************************************************************/
void Vorticity_2D_computation(Simulation *simu, int ifield_ux, int ifield_uy){
    field *f = simu->fd;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *temp_gradient;
    schnaps_real *vort;
    vort= calloc(nb_dof,sizeof(schnaps_real)); 
    temp_gradient= calloc(wdsize,sizeof(schnaps_real));
    //
    Compute_derivative(simu,temp_gradient,ifield_uy);
    //
    for (int i=0; i< nb_dof; i++){
      vort[i]= temp_gradient[i];
    };
    Compute_derivative(simu,temp_gradient,ifield_ux);
    for (int i=0; i< nb_dof; i++){
      vort[i]= vort[i]-temp_gradient[i+nb_dof];
    }; 
    //
    if (temp_gradient != NULL){
      free(temp_gradient);
    }
    if (vort != NULL){
      free(vort);
    }
}
/**************************************************************************************/
void Compute_and_dump_Vorticity_2d(Simulation *simu,char *filename, int create_file, schnaps_real t, int istep){
    LatticeData *ld=&schnaps_lattice_data;
    field *f = simu->fd;
    int ifield_ux=ld->index_ux;
    int ifield_uy=ld->index_uy;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *temp_gradient;
    schnaps_real *vort;
    vort= calloc(nb_dof,sizeof(schnaps_real)); 
    temp_gradient= calloc(wdsize,sizeof(schnaps_real));
    //
    Compute_derivative(simu,temp_gradient,ifield_uy);
    //
    for (int i=0; i< nb_dof; i++){
      vort[i]= temp_gradient[i];
    };
    Compute_derivative(simu,temp_gradient,ifield_ux);
    schnaps_real max_vort=-10000.0;
    for (int i=0; i< nb_dof; i++){
      vort[i]= vort[i]-temp_gradient[i+nb_dof];
      if (vort[i] > max_vort){
        max_vort=vort[i];
      }
    };
    printf("Max vort %f\n",max_vort); 
    //
    if (temp_gradient != NULL){
      free(temp_gradient);
    }
    //
    //PlotGenericFieldAsciiSparse(simu,vort,"vorticity",filename);
    //PlotGenericFieldBinSparse(simu,vort,"vorticity",filename);
    PlotExtScalarFieldBinMultitime(simu, vort,"vorticity",filename,create_file,t,istep);
    int typplot[3]= {ld->index_ux,ld->index_uy,ld->index_uz};
    PlotVecFieldsBinSparseMultitime(typplot,false,simu,"u",filename,0,t,istep);
    //
    if (vort != NULL){
      free(vort);
    }
}
/**************************************************************************************/
////// * test routine for field plot with sparser interpolation
//////////////////////////////////////
void PlotExtFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename) {
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
    simu->fd[0].raf[1],
    simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
    simu->fd[0].deg[1],
    simu->fd[0].deg[2]};
  const int npg[3] = {deg[0] + 1,
    deg[1] + 1,
    deg[2] + 1};
  const unsigned int sc_npg = npg[0] * npg[1] * npg[2];
  //
  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(schnaps_real));
  //
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  //int npgv = NPG(deg, nraf);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        }
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

        schnaps_real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};
        schnaps_real testpsi = 0;
        ////////////////////////////////////////
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projection to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
              schnaps_real psi;
              psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
              testpsi += psi;
              value[nodecount] += psi * w_in[index_glob_igp];
          }; // end loop subcell gauss points
        }; //end loop neighbour subcell 2
        };//end loop neighbour subcell 1
        }; //end loop neighbour subcell 0
        assert(fabs(testpsi-1) < _SMALL);
      ///////////////////////////////////////
      // Compare with an exact solution
      nodecount++;
      fprintf(gmshfile, "%d %f %f %f\n", nodecount,
        Xplot[0], Xplot[1], Xplot[2]);
      } // end loop Hex64 fe nodes
    } // end loop subcell 2
    } // end loop subcell 1
    } // end loop subcell 0
  } // end loop macrocells

  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  fprintf(gmshfile, "%d\n",
  simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]);
  int elm_type = 92;
  int num_tags = 0;
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
  	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    // Get the subcell id
    int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
    // Global subcell id
    int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
    fprintf(gmshfile, "%d ", numelem);
    fprintf(gmshfile, "%d ", elm_type);
    fprintf(gmshfile, "%d ", num_tags);
    for(int ii = 0; ii < 64; ii++) {
      int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
      fprintf(gmshfile, "%d ", numnoe);
    }
    fprintf(gmshfile, "\n");
    }
    }
    }
  } 
  fprintf(gmshfile, "$EndElements\n");
  // Now display data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field : anonymous\n");
  else
    fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);
  //
  for(int ino = 0; ino < nb_plotnodes; ino++) {
    fprintf(gmshfile, "%d %f\n", ino + 1, value[ino]);
  }
  fprintf(gmshfile, "\n$EndNodeData\n");
  //
  fclose(gmshfile);
  free(value);
}
void PlotExtFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename){
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  int one = 1;
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
            schnaps_real psi;
            psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
            testpsi += psi;
            value[nodecount] += psi * w_in[index_glob_igp];
          };
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        nodecount++;
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "field : anonymous");
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);

  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");

  fclose(gmshfile);
  free(value);
} 
/////////////////////////////////////
void PlotExtScalarFieldBinMultitime(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename,int create_file, schnaps_real t,int istep){
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  if (create_file==1){
  gmshfile = fopen(filename, "w" );
  }
  else
  {
  gmshfile = fopen(filename, "a" );
  }
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  if (create_file==1){
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  }
  int one = 1;
  if (create_file==1){
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  }
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  if (create_file==1){
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  }
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
            schnaps_real psi;
            psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
            testpsi += psi;
            value[nodecount] += psi * w_in[index_glob_igp];
          };
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        nodecount++;
        if (create_file==1){
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
        }
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  if (create_file==1){
  fprintf(gmshfile, "$EndNodes\n");
  }
  // Elements
  if (create_file==1){
    fprintf(gmshfile, "$Elements\n");
    int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
    fprintf(gmshfile, "%d\n",nb_elements);
    int elm_type = 92; //
    int num_tags = 0;
    int num_elm_follow= nb_elements;
    int elem_header[3] = {elm_type, num_elm_follow, num_tags};
    fwrite(elem_header, sizeof(int), 3, gmshfile);
    for(int i = 0; i < simu->macromesh.nbelems; i++) {
        // Loop on the subcells
        int icL[3];
        for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
        for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
        for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
        // Get the subcell id
        int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
        // Global subcell id
        int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
        int elm_data_size=1+num_tags+64; 
        int elem_data[elm_data_size];
        elem_data[0]= numelem;
          for(int ii = 0; ii < 64; ii++) {
            int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
            elem_data[ii+1]= numnoe;
          };
        fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
        };
        };
        };
    };
    fprintf(gmshfile, "$\nEndElements\n");
  }
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  int number_of_string_tags=1;
  fprintf(gmshfile, "%i\n",number_of_string_tags);
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field:anonymous \n");
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  int number_of_real_tags=1;
  fprintf(gmshfile,"%i\n",number_of_real_tags);
  fprintf(gmshfile,"%f\n",t);
  int number_of_int_tags=3;
  fprintf(gmshfile,"%i\n",number_of_int_tags);
  fprintf(gmshfile,"%i\n",istep);
  int number_of_components=1;
  fprintf(gmshfile,"%i\n",number_of_components);
  fprintf(gmshfile, "%i\n", nb_plotnodes);
  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");
  fclose(gmshfile);
  free(value);
} 
void PlotVecFieldsBinSparseMultitime(int typplot[3], int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep){

  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  if (create_file==1){
  gmshfile = fopen(filename, "w" );
  //printf("Creating file\n");
  }
  else
  {
  //printf("Appending data\n");
  gmshfile = fopen(filename, "a" );
  }
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  if (create_file==1){
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  }
  int one = 1;
  if (create_file==1){
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  }
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  if (create_file==1){
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  }
  //
  schnaps_real *value = calloc(3*nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
            schnaps_real psi;
            psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
            testpsi += psi;
            for (int idim=0;idim < 3; idim++){
            int vi = f->varindex(f->deg, f->raf, f->model.m, index_glob_igp, typplot[idim]);
            value[3*nodecount+idim] += psi * f->wn[vi];
            }
          };
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        // Compare with an exact solution
        if (compare) {
          schnaps_real wex[f->model.m];
          f->model.ImposedData(Xphy, f->tnow, wex);
          for (int idim=0;idim<3;idim++){
          value[3*nodecount+idim] -= wex[typplot[idim]];
          }
        }
        nodecount++;
        if (create_file==1){
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
        }
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  if (create_file==1){
  fprintf(gmshfile, "$EndNodes\n");
  }
  // Elements
  if (create_file==1){
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  }
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  int number_of_string_tags=1;
  fprintf(gmshfile, "%i\n",number_of_string_tags);
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", typplot);
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  int number_of_real_tags=1;
  fprintf(gmshfile,"%i\n",number_of_real_tags);
  fprintf(gmshfile,"%f\n",t);
  int number_of_int_tags=3;
  fprintf(gmshfile,"%i\n",number_of_int_tags);
  fprintf(gmshfile,"%i\n",istep);
  int number_of_components=3;
  fprintf(gmshfile,"%i\n",number_of_components);
  fprintf(gmshfile, "%i\n", nb_plotnodes);
  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    for (int idim=0; idim < 3;idim++){
    fwrite(&value[3*(ino-1)+idim],sizeof(double),1,gmshfile);
    }
  };
  fprintf(gmshfile, "\n$EndNodeData\n");
  fclose(gmshfile);
  free(value);
} 
//
void PlotFieldsBinSparseMultitime(int typplot, int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep){

  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  if (create_file==1){
  gmshfile = fopen(filename, "w" );
  //printf("Creating file\n");
  }
  else
  {
  //printf("Appending data\n");
  gmshfile = fopen(filename, "a" );
  }
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  if (create_file==1){
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  }
  int one = 1;
  if (create_file==1){
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  }
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  if (create_file==1){
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  }
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
            schnaps_real psi;
            psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
            testpsi += psi;
            int vi = f->varindex(f->deg, f->raf, f->model.m, index_glob_igp, typplot);
            value[nodecount] += psi * f->wn[vi];
          };
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        // Compare with an exact solution
        if (compare) {
          schnaps_real wex[f->model.m];
          f->model.ImposedData(Xphy, f->tnow, wex);
          value[nodecount] -= wex[typplot];
        }
        nodecount++;
        if (create_file==1){
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
        }
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  if (create_file==1){
  fprintf(gmshfile, "$EndNodes\n");
  }
  // Elements
  if (create_file==1){
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  }
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  int number_of_string_tags=1;
  fprintf(gmshfile, "%i\n",number_of_string_tags);
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", typplot);
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  int number_of_real_tags=1;
  fprintf(gmshfile,"%i\n",number_of_real_tags);
  fprintf(gmshfile,"%f\n",t);
  int number_of_int_tags=3;
  fprintf(gmshfile,"%i\n",number_of_int_tags);
  fprintf(gmshfile,"%i\n",istep);
  int number_of_components=1;
  fprintf(gmshfile,"%i\n",number_of_components);
  fprintf(gmshfile, "%i\n", nb_plotnodes);
  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");
  fclose(gmshfile);
  free(value);
} 
//**************************************************************************************************************************//
// lattice diagnostics (1D time traces) 
void Store_Lattice_diags(Simulation *simu){
  LatticeData *ld=&schnaps_lattice_data;
  schnaps_real tnow=simu->tnow;
  int iter= simu->iter_time_rk;
  int maxiter=simu->itermax_rk;
  int nb_diags=simu->nb_diags;
  //
  schnaps_real diag_vals[nb_diags];
  for (int irank=0;irank < nb_diags;irank++){
    diag_vals[irank]=0.0;
  };
  // actual collection
  ld->collect_diags(simu,diag_vals);
  //
  for (int irank=0;irank < nb_diags;irank++){
    simu->Diagnostics[iter* nb_diags+irank]= diag_vals[irank];
  };
}
void Dump_Lattice_Diagnostics(Simulation *simu,char simtag[3]){
  FILE *diagfile;
  assert((simu->Diagnostics != NULL));
  int nb_diags=simu->nb_diags;
  schnaps_real dt=simu->dt;
  char filename[sizeof("lbm_diag_TAG_raf000_cfl0.000.dat")];
  int raf=simu->fd[0].raf[0];
  schnaps_real cfl=simu->cfl;
  sprintf(filename,"lbm_diag_%s_raf%03d_cfl%1.3f.dat",simtag,raf,cfl);
  diagfile = fopen(filename,"w");
  for (int it=0; it < simu->itermax_rk;it++){
    schnaps_real tnow = dt * (schnaps_real) it;
    fprintf(diagfile,"%f\t",tnow);
    for (int irank=0;irank< simu->nb_diags;irank++){
      fprintf(diagfile,"%f\t",simu->Diagnostics[it*nb_diags+irank]);
    }
    fprintf(diagfile,"\n");
  };
  fclose(diagfile);
}
void Linear2DWave_CollectDiags(void *s,schnaps_real *diag_vals){
  Simulation *simu=s;
  schnaps_real error = 0;
  schnaps_real mean = 0;
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field *f = simu->fd + ie;
    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      schnaps_real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      w[iv] = f->wn[imem];
      }
      schnaps_real wex[f->model.m];
      schnaps_real wpg, det;
      // Compute wpg, det, and the exact solution
      schnaps_real xphy[3], xpgref[3];
      schnaps_real dtau[3][3], codtau[3][3];
      // Get the coordinates of the Gauss point
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
      xpgref, // xref
      NULL, -1, // dpsiref, ifa
      xphy, dtau, // xphy, dtau
      codtau, NULL, NULL); // codtau, dpsi, vnds
      det = dot_product(dtau[0], codtau[0]);
      // Get the exact value
      f->model.ImposedData(xphy, simu->tnow, wex);
    	schnaps_real diff = w[0] - wex[0];
      error += diff * diff * wpg * det;
      mean += wex[0] * wex[0] * wpg * det;
    };
    };
  //printf("errl2=%f\n",sqrt(error) / (sqrt(mean)  + 1e-14));
  diag_vals[0]=sqrt(error);
  diag_vals[1]=sqrt(mean);
  diag_vals[2]=sqrt(error) / (sqrt(mean)  + 1e-14);
}
void LatticeThetaTimeScheme(Simulation *simu, schnaps_real tmax, schnaps_real dt){
  LatticeData *ld=&schnaps_lattice_data;
  int nb_nodes=ld->q; // velocity nodes on the lattice
  int ig_glob=0,ig_node=0;
  field *f_glob,*f_node;
  // 
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};  
  int nb_ipg = NPG(deg, nraf);
  //
  Model model_advec;
  model_advec.m=1;
  model_advec.InitData =Lattice_Dummy_m1_InitData;
  model_advec.NumFlux=Lattice_OneNodeNumFlux;
  model_advec.ImposedData = Linear2DWave_ImposedData_OneNode;
  model_advec.BoundaryFlux = Linear2DWave_Periodic_BoundaryFlux_OneNode;
  model_advec.Source = NULL;
  Simulation *simu_advec;
  simu_advec=malloc(sizeof(Simulation));
  EmptySimulation(simu_advec);
  InitSimulation(simu_advec, &(simu->macromesh), deg, nraf, &model_advec);
  simu_advec->vmax=simu->vmax;
  simu_advec->cfl=0.0;
  simu_advec->nb_diags=0;
  simu_advec->pre_dtfields=NULL;
  simu_advec->post_dtfields=NULL;
  simu_advec->update_after_rk=NULL;
  //
  LinearSolver solver_imp[nb_nodes]; 
  LinearSolver solver_exp[nb_nodes];
  schnaps_real *res_advec=  calloc(simu_advec->wsize,sizeof(schnaps_real*)); 
  //
  schnaps_real theta=0.5; // implicit explicit splitting parameter 
  simu->dt=dt;
  int itermax=(int) (tmax/dt)+1;
  printf("Called with tmax=%f dt=%f Nb iterations:%i \n",tmax,dt,itermax);
  simu->itermax_rk=itermax;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  simu->tnow=0.0;
  simu_advec->dt=dt;
  simu_advec->itermax_rk=itermax;
  simu_advec->tnow=0.0;
  for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
    simu->fd[ie].tnow=simu->tnow;
    simu->fd[ie].tnow=simu_advec->tnow;
  }
  // Diagnostics (this should be elsewhere, some timetraces  module ?
  int size_diags = simu->nb_diags * itermax;
  if(simu->nb_diags != 0) {
     simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));
  };
  //
  time_t t_start,t_end; // time measurements for op factorization
  t_start=time(NULL);
  //  Solvers Init/Assembly
  printf("Sparse Linear Solvers init");
  for (int isim=0;isim < nb_nodes;isim++){
    ld->current_node_index=isim;
    //
    InitImplicitLinearSolver(simu_advec, &solver_imp[isim]);
    InitImplicitLinearSolver(simu_advec, &solver_exp[isim]);
    //
  }
  // End Operators init
  //
  // Time loop start
  for (int iter=0;iter < itermax;iter++){
    //
    if (iter%freq == 0){
      printf(" iter %i/%i t=%f\n",iter,itermax,simu->tnow);
    }
    //
    simu->iter_time_rk=iter;
    //
    if (iter==0){
      t_start=time(NULL);
    }
    //
    if (simu->pre_dtfields !=NULL){
      simu->pre_dtfields(simu);
    };
    // now loop on velocity nodes
    for(int isim=0;isim < nb_nodes;isim++){
      ld->current_node_index=isim;
      // dispatch main w to per node w's
      for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
        f_glob=simu->fd+ie;
        for(int ipg=0; ipg < nb_ipg;ipg++){
          f_node= simu_advec->fd+ie;
          ig_glob=f_glob->varindex(f_glob->deg,f_glob->raf,f_glob->model.m,ipg,isim);
          ig_node=f_node->varindex(f_node->deg,f_node->raf,1,ipg,0);
          f_node->wn[ig_node]= f_glob->wn[ig_glob];
        //
        }; // ipg end loop glops
      }; // ie end loop macroelements
      // end of data dispatch
      // Solvers assembly and factorization if necessary
      solver_imp[isim].rhs_is_assembly=false;
      solver_exp[isim].rhs_is_assembly=false;
      if (iter==0){
          solver_imp[isim].mat_is_assembly=false;
          solver_exp[isim].mat_is_assembly=false;
      }
      else
      {
          solver_imp[isim].mat_is_assembly=true;
          solver_exp[isim].mat_is_assembly=true;
      };
      //
      simu_advec->tnow=simu->tnow;
      for(int ie=0; ie < simu_advec->macromesh.nbelems; ++ie){
        simu_advec->fd[ie].tnow=simu_advec->tnow;
      }
      AssemblyImplicitLinearSolver(simu_advec, &solver_exp[isim],-(1.0-theta),simu_advec->dt);
      simu_advec->tnow=simu->tnow+simu->dt;
      for(int ie=0; ie < simu_advec->macromesh.nbelems; ++ie){
        simu_advec->fd[ie].tnow=simu_advec->tnow;
      }
      AssemblyImplicitLinearSolver(simu_advec, &solver_imp[isim],theta,simu_advec->dt);
      // compute residual
      MatVect(&solver_exp[isim], simu_advec->w, res_advec);
      //
      for(int i=0;i<solver_imp[isim].neq;i++){
        solver_imp[isim].rhs[i]=-solver_exp[isim].rhs[i]+solver_imp[isim].rhs[i]+res_advec[i];
      }
      //
      solver_imp[isim].solver_type=LU;
      Advanced_SolveLinearSolver(&solver_imp[isim],simu_advec);
      //
      for(int i=0;i<solver_imp[isim].neq;i++){
        simu_advec->w[i]=solver_imp[isim].sol[i];
      }
      // collect nodes ws to main w
      for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
        f_glob=simu->fd+ie;
        for(int ipg=0; ipg < nb_ipg;ipg++){
          f_node= simu_advec->fd+ie;
          ig_glob=f_glob->varindex(f_glob->deg,f_glob->raf,f_glob->model.m,ipg,isim);
          ig_node=f_node->varindex(f_node->deg,f_node->raf,1,ipg,0);
          f_glob->wn[ig_glob]= f_node->wn[ig_node]; 
        }; // ipg end loop glops
      }; // ie end loop macroelements
    }; // isim end loop on velocity node 
    // post advec ops
    if (simu->post_dtfields !=NULL){
      simu->post_dtfields(simu);
    }
    if (simu->update_after_rk !=NULL){
      simu->update_after_rk(simu,simu->w);
    }
    if (iter==0){
      t_end= time(NULL);
      printf("First step duration %f\n", t_end-t_start);
    }
    simu->tnow += simu->dt;
    simu_advec->tnow=simu->tnow;
    for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
      simu->fd[ie].tnow=simu->tnow;
      simu_advec->fd[ie].tnow=simu->tnow;
    }
    //
  }; // iter end of time loop 
  //
  if (res_advec !=NULL){
    free(res_advec);
  }
  if (simu_advec!= NULL){
    free(simu_advec);
  }
  //
}
/***************************************************************************************************************/