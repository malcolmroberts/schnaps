#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "model.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#define _XOPEN_SOURCE 700

double maxerr(double *a, double *b, int n) 
{
  double err = 0.0;
  for(int i = 0; i < n; ++i) {
    err = fmax(fabs(a[i] - b[i]), err);
  }
  return err;
}

int TestmEq2(void) {
  bool test = true;
  field f;
  
  f.model.cfl = 0.05;  
  f.model.m = 2; 
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  f.model.NumFlux = VecTransNumFlux2d;
  f.model.BoundaryFlux = VecTransBoundaryFlux2d;
  f.model.InitData = VecTransInitData2d;
  f.model.ImposedData = VecTransImposedData2d;
  f.varindex = GenericVarindex;

  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "VecTransNumFlux2d");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "VecTransBoundaryFlux2d");
  strcat(cl_buildoptions, buf);

  // Read the gmsh file
  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");

  // Try to detect a 2d mesh
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);

  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  Initfield(&f);
  f.is2d = true;
  
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  double *dtwn_cl = f.dtwn;
  double *dtwn = calloc(f.wsize, sizeof(double));
  
  double err;
  double tolerance = 1e-8;

  printf("Test volume terms\n");
  
  for(int i = 0; i < f.wsize; ++i) {
    dtwn_cl[i] = 0.0;
    dtwn[i] = 0.0;
  }
  
  f.dtwn = dtwn_cl;

  set_buf_to_zero_cl(&(f.dtwn_cl), f.wsize, &f);//, 0, NULL, NULL);
  clFinish(f.cli.commandqueue);

  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    update_physnode_cl(&f, ie, f.physnode_cl, f.physnode);
    //		       0, NULL, NULL);
    clFinish(f.cli.commandqueue);

    DGVolume_CL((void*) &(f.mcell[ie]), &f, &(f.wn_cl));//, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }

  CopyfieldtoCPU(&f);
  clFinish(f.cli.commandqueue);
  
  f.dtwn = dtwn;
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    DGSubCellInterface((void*) &(f.mcell[ie]), &f, f.wn, f.dtwn);
    DGVolume((void*) &(f.mcell[ie]), &f, f.wn, f.dtwn);
  }

  err = maxerr(dtwn, dtwn_cl, f.wsize);
  printf("\tmax error: %f\n", err);
  test = (err < tolerance);


  printf("Test macrocell interfaces\n");

  for(int i = 0; i < f.wsize; ++i) {
    dtwn_cl[i] = 0.0;
    dtwn[i] = 0.0;
  }

  f.dtwn = dtwn_cl;
  set_buf_to_zero_cl(&(f.dtwn_cl), f.wsize, &f);//, 0, NULL, NULL);
  clFinish(f.cli.commandqueue);
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ++ifa) {
    DGMacroCellInterface_CL((void*) (f.mface + ifa), &f, &f.wn_cl);
			    //			    0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }
  CopyfieldtoCPU(&f);
  clFinish(f.cli.commandqueue);
  
  f.dtwn = dtwn;
  MacroFace mface[f.macromesh.nbfaces];
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    mface[ifa].first = ifa;
    mface[ifa].last_p1 = ifa + 1;
  }
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    DGMacroCellInterface((void*) (mface + ifa), &f, f.wn, f.dtwn);
  }
  err = maxerr(dtwn, dtwn_cl, f.wsize);
  printf("\tmax error: %f\n", err);
  test = (err < tolerance);

  printf("Test all terms\n");
  for(int i = 0; i < f.wsize; ++i) {
    dtwn_cl[i] = 0.0;
    dtwn[i] = 0.0;
  }
  
  f.dtwn = dtwn_cl;
  dtfield_CL(&f, &(f.wn_cl));//, 0, NULL, NULL);
  clFinish(f.cli.commandqueue);
  
  CopyfieldtoCPU(&f);
  clFinish(f.cli.commandqueue);
  
  f.dtwn = dtwn;
  dtfield(&f, f.wn, f.dtwn);

  err = maxerr(dtwn, dtwn_cl, f.wsize);
  printf("\tmax error: %f\n", err);
  test = (err < tolerance);

  //Displayfield(&f);
 
  /* double tmax = 0.1; */
  /* //RK2_CL(&f, tmax, 0, NULL, NULL); */
  /* CopyfieldtoCPU(&f); */
 
  /* // Save the results and the error */
  /* Plotfield(0, false, &f, NULL, "dgvisu.msh"); */
  /* Plotfield(0, true, &f, "error", "dgerror.msh"); */

  /* double dd = L2error(&f); */
  /* double tolerance = 1e-4; */
  /* test = test && (dd < tolerance); */
  /* printf("L2 error: %f\n", dd); */

  return test;
};

int main(void) {
  int resu = TestmEq2();
  if (resu) 
    printf("OpenCL m greater than 1 test OK!\n");
  else 
    printf("OpenCL m greater than 1 test failed !\n");
  return !resu;
}
