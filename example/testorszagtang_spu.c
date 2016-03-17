#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "../test/test.h"
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
int TestOrszagTang_SPU(int argc, char *argv[]) ;


int main(int argc, char *argv[]) {
  int resu = TestOrszagTang_SPU(argc,argv);
  if (resu)
    printf("OrszagTang test OK !\n");
  else 
    printf("OrszagTang test failed !\n");
  return !resu;
}

int TestOrszagTang_SPU(int argc, char *argv[]) {

  int test = true;

  #ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
  #endif


  // disable OPENCL for StarPU
  putenv("STARPU_NOPENCL=0");

  ////////////////////
  // init mesh
  MacroMesh mesh;
  char *mshname =  "../test/testOTgrid.msh";
  
  ReadMacroMesh(&mesh, mshname);
  Detect2DMacroMesh(&mesh);
  bool is2d=mesh.is2d; 
  assert(is2d);  

  schnaps_real periodsize = 6.2831853;
  mesh.period[0]=periodsize;
  mesh.period[1]=periodsize;

  BuildConnectivity(&mesh);

  int deg[]={2, 2, 0};
  int raf[]={15, 15, 1};
  CheckMacroMesh(&mesh, deg, raf);

  ////////////////////

  
  Model model;

  Simulation simu;
  EmptySimulation(&simu);

  schnaps_real cfl = 0.2;
  simu.cfl = cfl;
  model.m = 9;

  strcpy(model.name,"MHD");

  model.NumFlux=MHDNumFluxRusanov;
  model.BoundaryFlux=MHDBoundaryFluxOrszagTang;
  model.InitData=MHDInitDataOrszagTang;
  model.ImposedData=MHDImposedDataOrszagTang;
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

  sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math", "MHDBoundaryFluxOrszagTang");
  strcat(cl_buildoptions, buf);


  starpu_use = true;
  starpu_c_use = true;
  starpu_ocl_use = true;
  
  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  schnaps_real tmax = 1.0;
  simu.vmax = 6.0;
  schnaps_real dt = 0;
  
  RK2_SPU(&simu, tmax);
  //RK4(&simu, tmax);
  
  //CopyfieldtoCPU(&simu);
 
  //show_cl_timing(&simu);
  PlotFields(0, false, &simu, NULL, "dgvisu.msh");

  return test;
}
