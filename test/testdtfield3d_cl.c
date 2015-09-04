#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

int TestDtfield3D_CL(void)
{
  return false; // FIXME!!!!
  
  /* bool test = true; */

  /* if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) { */
  /*   printf("OpenCL device not acceptable.\n"); */
  /*   return true; */
  /* } */

  /* field f; */
  /* init_empty_field(&f); */


  /* f.model.cfl = 0.05;   */
  /* f.model.m = 8; // num of conservative variables */

  /* f.model.NumFlux = Maxwell3DNumFluxClean_uncentered; */
  /* f.model.BoundaryFlux = Maxwell3DBoundaryFlux_uncentered; */
  /* f.model.InitData = Maxwell3DInitData; */
  /* f.model.ImposedData = Maxwell3DImposedData; */
  /* f.varindex = GenericVarindex; */
    
  /* f.interp.interp_param[0] = f.model.m; */
  /* f.interp.interp_param[1] = 3; // x direction degree */
  /* f.interp.interp_param[2] = 3; // y direction degree */
  /* f.interp.interp_param[3] = 3; // z direction degree */
  /* f.interp.interp_param[4] = 4; // x direction refinement */
  /* f.interp.interp_param[5] = 4; // y direction refinement */
  /* f.interp.interp_param[6] = 4; // z direction refinement */

  /* ReadMacroMesh(&(f.macromesh), "../test/testcube.msh"); */

  /* // FIXME: temp */
  /* /\* Detect2DMacroMesh(&(f.macromesh)); *\/ */
  /* /\* assert(f.macromesh.is2d); *\/ */

  /* BuildConnectivity(&(f.macromesh)); */

  /* char buf[1000]; */
  /* sprintf(buf, "-D _M=%d", f.model.m); */
  /* strcat(cl_buildoptions, buf); */

  /* // Source is for rho and J, which are zero here. */
  /* //set_source_CL(&f, "Maxwell3DSource"); */
  /* sprintf(numflux_cl_name, "%s", "Maxwell3DNumFluxClean_uncentered"); */
  /* sprintf(buf," -D NUMFLUX="); */
  /* strcat(buf, numflux_cl_name); */
  /* strcat(cl_buildoptions, buf); */

  /* sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell3DBoundaryFlux_uncentered"); */
  /* strcat(cl_buildoptions, buf); */

  /* Initfield(&f); */
  
  /* cl_event clv_dtfield = clCreateUserEvent(f.cli.context, NULL); */
  
  /* dtfield_CL(&f, &f.wn_cl, 0, NULL, &clv_dtfield); */
  /* clWaitForEvents(1, &clv_dtfield); */
  /* CopyfieldtoCPU(&f); */

  /* // Displayfield(&f); */
  /* show_cl_timing(&f); */

  /* real *saveptr = f.dtwn; */
  /* f.dtwn = calloc(f.wsize, sizeof(real)); */

  /* //f.model.Source = OneSource; */
  /* dtfield(&f, f.wn, f.dtwn); */
 
  /* real maxerr = 0; */
  /* for(int i = 0; i < f.wsize; i++) { */
  /*   real error = f.dtwn[i] - saveptr[i]; */
  /*   //printf("error= \t%f\t%f\t%f\n", error, f.dtwn[i], saveptr[i]); */
  /*   maxerr = fmax(fabs(error), maxerr); */
  /* } */
  /* printf("max error: %f\n", maxerr); */

  /* test = (maxerr < 0.00005); */

  /* return test; */
}

int main(void) {
  int resu = TestDtfield3D_CL();

  if (resu) 
    printf("3D Dtfield_CL test OK !\n");
  else 
    printf("3D Dtfield_CL test failed !\n");

  return !resu;
} 
