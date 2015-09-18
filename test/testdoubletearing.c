#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "mhd.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#define _XOPEN_SOURCE 700

//real seconds()
//{
//  struct timespec ts;
//  clock_gettime(CLOCK_MONOTONIC, &ts);
//  return (real)ts.tv_sec + 1e-9 * (real)ts.tv_nsec;
//}

int main(int argc, char *argv[]) {
  int resu = TestDoubleTearing(argc,argv);
  if (resu)
    printf("DoubleTearing test OK !\n");
  else 
    printf("DoubleTearing test failed !\n");
  return !resu;
}

int TestDoubleTearing(int argc, char *argv[])
{

  int test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  Simulation simu;
  EmptySimulation(&simu);

  MacroMesh mesh;
  char *mshname =  "../test/testdoubletearinggrid.msh";

  ReadMacroMesh(&mesh, mshname);
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);
  int deg[]={1, 1, 0};
  int raf[]={10, 10, 1};
  CheckMacroMesh(&mesh, deg, raf);

  Model model;

  assert(mesh.is2d);

  simu.cfl = 0.1;
  model.m = 9;
  real periodsize = 4.;

  model.NumFlux = MHDNumFluxRusanov;
  model.BoundaryFlux = MHDBoundaryFluxDoubleTearing;
  model.InitData = MHDInitDataDoubleTearing;
  model.ImposedData = MHDImposedDataDoubleTearing;
  model.Source = NULL;
  
  char buf[1000];
  sprintf(buf, "-D _M=%d -D _PERIODX=%f -D _PERIODY=%f",
          model.m,
          periodsize,
          periodsize);
  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", "MHDNumFluxRusanov");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "MHDBoundaryFluxDoubleTearing");
  strcat(cl_buildoptions, buf);

  simu.macromesh.period[1]=periodsize;
  
  //AffineMapMacroMesh(&(simu.macromesh));
  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  real tmax = 0.1;
  simu.vmax = 6.0;
  real dt = 0;
  RK2_CL(&simu, tmax, dt,  0, NULL, NULL);
  
  CopyfieldtoCPU(&simu);
 
  PlotFields(2, false, &simu, "P", "dgvisu.msh");
  //Plotfield(0, true , &simu, "error", "dgerror.msh");

  real dd = L2error(&simu);

  printf("L2 error: %f\n", dd);

  show_cl_timing(&simu);

  return test;

}
