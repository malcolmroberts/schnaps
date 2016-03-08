#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

int TestfieldRK2_SPU(void){
  bool test = true;

  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh


  putenv("STARPU_NOPENCL=0");


  char *mshname =  "../test/disque2d.msh";

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testdisque.msh");
  //ReadMacroMesh(&mesh,"../test/testcube2.msh");
  //Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

  Simulation simu;
  EmptySimulation(&simu);

#if 0
  // 2D version
  model.cfl = 0.05;
  model.m = 1;

  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
  model.Source = NULL;

  char buf[1000];
#ifdef _DOUBLE_PRECISION
  sprintf(buf, "-D real=double -D _M=%d", model.m);
#else
  sprintf(buf, "-D real=float -D _M=%d", model.m);
#endif
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "TransNumFlux2d");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "TransBoundaryFlux2d");
  strcat(cl_buildoptions, buf);

  int deg[]={3, 3, 0};
  int raf[]={3, 3, 1};

  assert(mesh.is2d);
  assert(1==2);
#else
  // 3D version
  model.m = 1;
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;
  model.Source = NULL;

  char buf[1000];
#ifdef _DOUBLE_PRECISION
  sprintf(buf, "-D real=double -D _M=%d", model.m);
#else
  sprintf(buf, "-D real=float -D _M=%d", model.m);
#endif
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "TransNumFlux");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "TestTransBoundaryFlux");
  strcat(cl_buildoptions, buf);

  int deg[]={2, 2, 2};
  int raf[]={2, 2, 2};
  //int raf[]={2, 2, 2};

#endif

  // 2015-01-19: the below parameters fail with testmacrocellinterface
  // but pass the test here (perhaps because the error is hidden by
  // the RK error?)
  /*
  f.model.cfl = 0.05;
  f.model.m = 1;
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 1; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement
  */

  //AffineMapMacroMesh(&(f.macromesh));

#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif

  CheckMacroMesh(&mesh, deg, raf);
  starpu_use = true;
  starpu_c_use = true;
  starpu_ocl_use = true;

  InitSimulation(&simu, &mesh, deg, raf, &model);

  schnaps_real tmax = 0.25;
  tmax = 0.00666;
  simu.cfl=0.2;
  simu.vmax=1;
  RK2_SPU(&simu,tmax);

  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%.12f\n", dd);

  schnaps_real tolerance = 0.0026;

  test = dd < tolerance;

  FreeMacroMesh(&mesh);

  return test;
}

int main(void) {
  int resu = TestfieldRK2_SPU();
  if(resu)
    printf("starpu field RK2 test OK !\n");
  else
    printf("starpu field RK2 test failed !\n");
  return !resu;
}
