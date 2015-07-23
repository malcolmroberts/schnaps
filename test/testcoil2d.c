#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "maxwell.h"

void Coil2DImposedData(const real x[3],const real t,real w[])
{
  real r = x[0] * x[0] + x[1] * x[1];
  w[0] = 0;
  w[1] = 0;
  w[2] = r > 1 ? 0 : 1;
  w[3] = 0;
  w[4] = 0;
  w[5] = 0;
  w[6] = 0;
}

void coil_pre_dtfields(void *simu, real *w);

void coil_pre_dtfields(void *simu, real *w){
  AccumulateParticles(simu, w);
}


void Coil2DSource(const real *x, const real t, const real *w, real *source)
{
  // w: (Ex, Ey, Hz, Hz, \lambda, rho, Jx, Jy)
  
  // FIXME add documentation
  
  const real khi = 1.0;
  source[0] = -w[4];
  source[1] = -w[5];
  source[2] = 0;
  source[3] = 0;//instead of khi * w[6]: we want div E =0
  source[4] = 0;
  source[5] = 0;
  source[6] = 0;
}




void Coil2DBoundaryFlux(real x[3], real t, real wL[], real *vnorm,
			real *flux)
{
  real wR[7];
  Coil2DImposedData(x, t, wR);
  Maxwell2DNumFlux_uncentered(wL, wR, vnorm, flux);
}

void Coil2DInitData(real x[3], real w[])
{
  real t = 0;
  Coil2DImposedData(x, t, w);
}

int TestCoil2D(void)
{
  bool test = true;

  char *mshname =  "test/testmacromesh.msh";
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,mshname);
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);


  // test gmsh file reading
  ReadMacroMesh(&mesh, "test/testmacromesh.msh");
  Detect2DMacroMesh(&mesh);
  assert(mesh.is2d);
  BuildConnectivity(&mesh);
  //PrintMacroMesh(&m);

  Model model;

  model.m = 7; // num of conservative variables

  model.NumFlux = Maxwell2DNumFlux_uncentered;
  model.BoundaryFlux = Coil2DBoundaryFlux;
  model.InitData = Coil2DInitData;
  model.Source = Coil2DSource;
  model.ImposedData = Coil2DImposedData;
    
  int deg[]={2, 2, 0};
  int raf[]={8, 8, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.pre_dtfields = coil_pre_dtfields; // must be called after init

   
  
  PIC pic;
  simu.pic = &pic;

  InitPIC(&pic,100); 
  CreateCoil2DParticles(&pic, &mesh);
  PlotParticles(&pic, &mesh);

  // time evolution
  real tmax = 0.1;
  simu.cfl=0.2;
  simu.vmax = 1;
  RK2(&simu, tmax);
 
  // Save the results and the error
  PlotFields(2, false, &simu, NULL, "dgvisu.msh");
  PlotFields(2, true, &simu, "error", "dgerror.msh");

  real dd = L2error(&simu);
  real tolerance = 0.3;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
}

int main(void) {
  int resu = TestCoil2D();
  if (resu) 
    printf("Coil2D test OK!\n");
  else 
    printf("Coil2D failed !\n");
  return !resu;
}
