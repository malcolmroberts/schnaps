#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "clutils.h"

int TestMacroFace(void);

int main(void) {
  int resu = TestMacroFace();
  if(resu) 
    printf("MacroFace test OK\n");
  else 
    printf("MacroFace test FAILED\n");
  return !resu;
} 

int TestMacroFace(void){
  bool test = true;

  Field f;
  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 1; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testdisque.msh");
  BuildConnectivity(&(f.macromesh));
  InitField(&f);
  

  // OpenCL method
  // NB: InitField expects a certain address for dtwn, so the OpenCL
  // version must come before the other versions.
  cl_int status;
  void* chkptr = clEnqueueMapBuffer(f.cli.commandqueue,
				    f.dtwn_cl,
				    CL_TRUE,
				    CL_MAP_WRITE,
				    0, // offset
				    sizeof(cl_double) * 60,
				    0, NULL, NULL,
				    &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  assert(chkptr == f.dtwn);

  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;

  status = clEnqueueUnmapMemObject(f.cli.commandqueue,
				   f.dtwn_cl,
				   f.dtwn,
				   0, NULL, NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  clFinish(f.cli.commandqueue);  // wait the end of the computation

  CopyFieldtoCPU(&f);
  double *fdtwn_opencl = f.dtwn;

  // OpenMP, new method
  MacroFace mface[f.macromesh.nbfaces];
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    mface[ifa].field = &f;
    mface[ifa].first = ifa;
    mface[ifa].last_p1 = ifa + 1;
  }
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++)
    DGMacroCellInterface_CL((void*) (mface + ifa));
  f.dtwn = calloc(f.wsize, sizeof(double));
  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++)
    DGMacroCellInterface((void*) (mface + ifa));
  double *fdtwn_openmp = f.dtwn;

  // OpenMP, slow method
  MacroCell mcell[f.macromesh.nbelems];
  for(int ie = 0; ie < f.macromesh.nbelems; ie++) {
    mcell[ie].field = &f;
    mcell[ie].first = ie;
    mcell[ie].last_p1 = ie + 1;
  }
  f.dtwn = calloc(f.wsize, sizeof(double));
  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    MacroCell *mcelli = mcell + ie;
    DGMacroCellInterfaceSlow(mcelli);
  }
  double *fdtwn_slow = f.dtwn;

  // Check that the results are the same
  test = true;
  double tolerance = 1e-8;

  double maxerr = 0.0;
  for(int i = 0; i < f.wsize; i++)
    maxerr = fmax(fabs(fdtwn_openmp[i] - fdtwn_opencl[i]), maxerr);
  printf("Max difference between OpenCL and OpenMP: %f\n", maxerr);
  test = test && (maxerr < tolerance);

  maxerr = 0.0;
  for(int i = 0; i < f.wsize; i++)
    maxerr = fmax(fabs(fdtwn_openmp[i] - fdtwn_slow[i]), maxerr);
  printf("Max difference between OpenMP and OpenMP-slow: %f\n", maxerr);
  test = test && (maxerr < tolerance);

  maxerr = 0.0;
  for(int i = 0; i < f.wsize; i++)
    maxerr = fmax(fabs(fdtwn_opencl[i] - fdtwn_slow[i]), maxerr);
  printf("Max difference between OpenCL and OpenMP-slow: %f\n", maxerr);
  test = test && (maxerr < tolerance);

  return test;
}
