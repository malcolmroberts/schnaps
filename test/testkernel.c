#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestKernel(void)
{
  bool test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  Model model;
  model.m = m;
  
  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
  
  /* field f; */
  /* init_empty_field(&f); */

  /* f.varindex = GenericVarindex; */

  int deg[3] = {1, 1, 0}; // Poynomial degree
  int raf[3] = {2, 2, 1}; // Number of subcells per macrocell
  
  /* f.model.cfl = 0.05; */


  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../test/testmacromesh.msh");
  Detect2DMacroMesh(&mesh);
  assert(mesh.is2d);
  BuildConnectivity(&mesh);

  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);

   // Set dtwn to 1 for testing
  {
    void *chkptr;
    cl_int status;
    chkptr = clEnqueueMapBuffer(simu.cli.commandqueue,
  				simu.dtw_cl,  // buffer to copy from
  				CL_TRUE,  // block until the buffer is available
  				CL_MAP_WRITE,
  				0, // offset
  				sizeof(real) * (simu.wsize), // buffersize
  				0,
  				NULL,
  				NULL, // events management
  				&status);
    assert(status == CL_SUCCESS);
    assert(chkptr == simu.dtw);

    for(int i = 0; i < simu.wsize; i++) {
      simu.dtw[i] = 1;
    }

    status = clEnqueueUnmapMemObject(simu.cli.commandqueue,
  				     simu.dtw_cl,
  				     simu.dtw,
  				     0,
  				     NULL,
  				     NULL);
    assert(status == CL_SUCCESS);

    status = clFinish(simu.cli.commandqueue);
    assert(status == CL_SUCCESS);
  }
 
  for(int ie = 0; ie < mesh.nbelems; ++ie) {
    printf("ie: %d\n", ie);
    DGMass_CL(ie, &simu, 0, NULL, NULL);
    clFinish(simu.cli.commandqueue);
  }

  CopyfieldtoCPU(&simu);

  /* Displayfield(&f); */

  // save the dtwn pointer
  real *saveptr = simu.dtw;

  // malloc a new dtwn.
  simu.dtw = calloc(simu.wsize, sizeof(real));
  for(int i = 0; i < simu.wsize; i++){
    simu.dtw[i] = 1;
  }

  int fsize =  simu.wsize / mesh.nbelems;
  for(int ie = 0; ie < mesh.nbelems; ++ie)
    DGMass(simu.fd + ie, simu.w + ie * fsize, simu.dtw + ie * fsize);

  //Displayfield(&f);

  //check that the results are the same
  real maxerr = 0;
  for(int i = 0; i < simu.wsize; i++){
    //printf("error=%f %f %f\n", f.dtw[i] - saveptr[i], f.dtw[i], saveptr[i]);
    maxerr = fmax(fabs(simu.dtw[i] - saveptr[i]), maxerr);
  }
  printf("max error=%f\n",maxerr);

  test = (maxerr < _SMALL);

  return test;
}

int main() {
  // Unit tests
  int resu = TestKernel();
  if (resu) 
    printf("Kernel test OK !\n");
  else 
    printf("Kernel test failed !\n");
  return !resu;
} 
