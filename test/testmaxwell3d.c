#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include <string.h>

int TestMaxwell3D(void) {
  bool test = true;

  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh

  char *mshname =  "../test/testcube.msh";
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  //ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  //Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

  model.m = 8;

  model.NumFlux = Maxwell3DNumFluxClean_upwind;
  //f.model.NumFlux = Maxwell2DNumFlux_centered;
  model.BoundaryFlux = Maxwell3DBoundaryFlux_upwind;
  model.InitData = Maxwell3DInitData;
  model.ImposedData = Maxwell3DImposedData;
  //model.Source = Maxwell2DSource;
  model.Source = NULL;


  int deg[]={2, 2, 2};
  int raf[]={32, 32, 32};

  //assert(mesh.is2d);

#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif
  
  CheckMacroMesh(&mesh, deg, raf);



  Simulation simu;
  EmptySimulation(&simu);

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  //set_source_CL(&simu, "Maxwell3DSource");
  sprintf(numflux_cl_name, "%s", "Maxwell3DNumFluxClean_upwind");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell3DBoundaryFlux_upwind");
  strcat(cl_buildoptions, buf);



  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  real tmax = .1;
  simu.cfl=0.2;
  simu.vmax=1;

#if 0
  // C version
  RK4(&simu, tmax);
#else
  // OpenCL version
  real dt = 0;
  RK4_CL(&simu, tmax, dt, 0, 0, 0);

  CopyfieldtoCPU(&simu); 
  printf("\nOpenCL Kernel time:\n");
  show_cl_timing(&simu);
  printf("\n");
#endif


  //PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  //PlotFields(0, true , &simu, "error", "dgerror.msh");

  real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%f\n", dd);

  real tolerance = 0.0025;
  tolerance = 0.08;
  test = dd < tolerance;
  
 
  return test;
}

int main(void) {
  int resu = TestMaxwell3D();
  if (resu) 
    printf("Maxwell3D test OK!\n");
  else 
    printf("Maxwell3D failed !\n");
  return !resu;
}
