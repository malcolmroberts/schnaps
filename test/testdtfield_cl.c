#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

int TestDtfield_CL(void);

int main(void) {
  int resu = TestDtfield_CL();

  if (resu) 
    printf("Dtfield_CL test OK !\n");
  else 
    printf("Dtfield_CL test failed !\n");

  return !resu;
} 

int TestDtfield_CL(void){
  bool test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  
  MacroMesh mesh;
  char *mshname = "../test/disque2d.msh";
  //char *mshname = "../test/unit-cube.msh";
  
  ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  //ReadMacroMesh(&mesh,"../test/testcube2.msh");
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

#if 1
  // 2D version
  assert(mesh.is2d);

  model.cfl = 0.05;
  model.m = 1;
  m = model.m;

  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
 
  int deg[]={2, 2, 0};
  int raf[]={4, 4, 1};

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "TransNumFlux2d");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "TransBoundaryFlux2d");
  strcat(cl_buildoptions, buf);


#else
  // 3D version
  assert(!mesh.is2d);
  
  model.cfl = 0.05;
  model.m = 1;
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;

  int deg[]={2, 2, 2};
  int raf[]={3, 3, 3};
  
#endif

  set_global_m(model.m);

  Simulation simu;

  model.Source = OneSource;
  set_source_CL(&simu, "OneSource");
  //model.Source = NULL;


  assert(simu.use_source_cl);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  cl_int status;
  cl_event clv_dtfield;// = clCreateUserEvent(simu.cli.context, &status);
  printf("dtfield_CL...\n");
  dtfield_CL(&simu, &simu.w_cl, 0, NULL, &clv_dtfield);
  printf("\t...dtfield_CL done\n");
  
  //clWaitForEvents(1, &clv_dtfield);
  CopyfieldtoCPU(&simu);

  // Displayfield(&f);
  show_cl_timing(&simu);

  real *saveptr = simu.dtw;
  simu.dtw = calloc(simu.wsize, sizeof(real));

  DtFields(&simu, simu.w, simu.dtw);
 
  real maxerr = 0;
  for(int i = 0; i < simu.wsize; i++) {
    real error = simu.dtw[i] - saveptr[i];
    //printf("error= \t%f\t%f\t%f\n", error, f.dtwn[i], saveptr[i]);
    maxerr = fmax(fabs(error), maxerr);
  }
  printf("max error: %f\n", maxerr);

  test = (maxerr < _SMALL * 10);

  return test;
}
