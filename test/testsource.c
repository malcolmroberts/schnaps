#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "clutils.h"

void testSource(const real *x, const real t, const real *w, real *source,
		int m) 
{
  for(int i = 0; i < m; ++i) {
    source[i] = 1.0;
  }
}

int TestKernel()
{
  int retval = 0;

  real tolerance;
  if(sizeof(real) == sizeof(double))
    tolerance = 1e-8;
  else
    tolerance = 1e-4;
  
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  field f;
  init_empty_field(&f);

  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;
  

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 1; // x direction degree
  f.interp.interp_param[2] = 1; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&f.macromesh, "../test/testmacromesh.msh");
  //ReadMacroMesh(&f.macromesh,"test/testcube.msh");
  Detect2DMacroMesh(&f.macromesh);
  assert(f.macromesh.is2d);  
  BuildConnectivity(&f.macromesh);

  set_source_CL(&f, "OneSource"); // MUST be done before Initfield
  Initfield(&f);
  f.model.Source = testSource; // MUST be done after Initfield

    
  for(int i = 0; i < f.wsize; i++){
    f.dtwn[i] = 0;
  }

  CopyfieldtoGPU(&f);

  real tnow = 0;
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    MacroCell *mcell = f.mcell + ie;
    DGSource_CL(mcell, &f, tnow, f.wn_cl + ie, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
    DGMass_CL(mcell, &f, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }

  CopyfieldtoCPU(&f);

  real maxerr = 0;
  for(int i = 0; i < f.wsize; i++) {
    real err = fabs(f.dtwn[i] - 1.0);
    maxerr = fmax(err, maxerr);
  }
  printf("\nOpenCL max error:\t%f\n", maxerr);

  if(maxerr > tolerance)
    retval += 1;
  
  //Displayfield(&f);

  // save the dtwn pointer
  real *saveptr = f.dtwn;

  // malloc a new dtwn.
  f.dtwn = calloc(f.wsize, sizeof(real));
  for(int i = 0; i < f.wsize; i++){
    f.dtwn[i] = 0;
  }

  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    MacroCell *mcell = f.mcell + ie;
    real *dtwnmc = f.dtwn + mcell->woffset;
    real *wmc = f.wn + mcell->woffset;
    DGSource(mcell, &f, tnow, wmc, dtwnmc);
    DGMass(mcell, &f, dtwnmc);
  }
  
  maxerr = 0;
  for(int i = 0; i < f.wsize; i++) {
    real err = fabs(f.dtwn[i] - 1.0);
    maxerr = fmax(err, maxerr);
  }
  printf("\nC max error:      \t%f\n",maxerr);

  if(maxerr > tolerance)
    retval += 1;

  assert(f.dtwn != saveptr);
  maxerr = 0;
  for(int i = 0; i < f.wsize; i++){
    real err = fabs(f.dtwn[i] - saveptr[i]);
    maxerr = fmax(err, maxerr);
  }
  printf("\nmax difference:      \t%f\n",maxerr);

  if(maxerr > tolerance)
    retval += 1;
  
  return retval;
}

int main()
{
  int retval = TestKernel();
  if (retval == 0) 
    printf("Kernel test OK !\n");
  else 
    printf("Kernel test failed !\n");
  return retval;
} 
