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
  int resu = TestReconnexion(argc,argv);
  if (resu)
    printf("Reconnexion test OK !\n");
  else 
    printf("Reconnexion test failed !\n");
  return !resu;
}

int TestReconnexion(int argc, char *argv[])
{
  
  int test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  Simulation simu;
  EmptySimulation(&simu);

  MacroMesh mesh;
  char *mshname =  "../test/testrecogrid.msh";

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
  real periodsize = 2.;

  model.NumFlux = MHDNumFluxRusanov;
  model.BoundaryFlux = MHDBoundaryFluxReconnexion;
  model.InitData = MHDInitDataReconnexion;
  model.ImposedData = MHDImposedDataReconnexion;
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

  sprintf(buf, " -D BOUNDARYFLUX=%s", "MHDBoundaryFluxReconnexion");
  strcat(cl_buildoptions, buf);

  simu.macromesh.period[0]=periodsize;
  simu.macromesh.period[1]=periodsize;
  
  //AffineMapMacroMesh(&(simu.macromesh));
  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  real tmax = 0.1;
  simu.vmax = 6.0;
  real dt = 0;
  RK2_CL(&simu, tmax, dt,  0, NULL, NULL);
  
  CopyfieldtoCPU(&simu);
 
  PlotFields(5, false, &simu, "By", "dgvisu.msh");
  //Plotfield(0, true , &simu, "error", "dgerror.msh");

  real dd = L2error(&simu);

  printf("L2 error: %f\n", dd);

  show_cl_timing(&simu);

  return test;

}
