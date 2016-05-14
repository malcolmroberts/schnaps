#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include <string.h>

int TestMaxwell2D()
{
  bool test = true;
  field f;
  init_empty_field(&f);

  f.model.cfl = 0.05;  
  f.model.m = 7; // num of conservative variables

  f.model.NumFlux = Maxwell2DNumFlux_upwind;
  //f.model.NumFlux = Maxwell2DNumFlux_centered;
  f.model.BoundaryFlux = Maxwell2DBoundaryFlux_upwind;
  f.model.InitData = Maxwell2DInitData;
  f.model.ImposedData = Maxwell2DImposedData;
  f.varindex = GenericVarindex;
  f.model.Source = Maxwell2DSource;
  
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "../test/testcube.msh");

  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);

  BuildConnectivity(&(f.macromesh));

  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  set_source_CL(&f, "Maxwell2DSource");
  sprintf(numflux_cl_name, "%s", "Maxwell2DNumFlux_upwind");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell2DBoundaryFlux_upwind");
  strcat(cl_buildoptions, buf);

  Initfield(&f);
  
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  real tmax = 0.1;
  f.vmax = 1;
  real dt = set_dt(&f);

#if 0
  // C version
  RK2(&f, tmax, dt);
#else
  // OpenCL version
  CopyfieldtoGPU(&f);
  RK2_CL(&f, tmax, dt, 0, 0, 0);
  CopyfieldtoCPU(&f);
  printf("\nOpenCL Kernel time:\n");
  show_cl_timing(&f);
  printf("\n");
#endif

  // Save the results and the error
  /* Plotfield(0, false, &f, NULL, "dgvisu.msh"); */
  /* Plotfield(0, true, &f, "error", "dgerror.msh"); */

  real dd = L2error(&f);
  real tolerance = 1.1e-2;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
}

int main(void) {
  int resu = TestMaxwell2D();
  if (resu) 
    printf("Maxwell2D test OK!\n");
  else 
    printf("Maxwell2D failed !\n");
  return !resu;
}
