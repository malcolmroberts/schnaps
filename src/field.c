#include "field.h"
#include "geometry.h"
#include "interpolation.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "global.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include "quantities_vp.h"
#include "solverpoisson.h"
#include "model.h"

#ifdef _WITH_OPENCL 
#include "clutils.h"
#include "clinfo.h"
#endif

#pragma start_opencl
int GenericVarindex(__constant int *deg, __constant int *raf,
		    const int m, int elem,
		    int ipg, int iv) {
  int npg = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1)
    * raf[0] * raf[1] * raf[2];
  return iv + m * (ipg + npg * elem);
}
#pragma end_opencl

#pragma start_opencl
int GenericVarindex3d(int m, int *nx, int *nc,
		      int elem,
		      int iv, int *ix, int *ic) 
{
  // NB: passing nx and nx separately may allow one to deal with faces better.
  // Also, keeping the values in registers is theoretically faster than
  // even __local memory.

  // int nx[3] = {param[1] + 1, param[2] + 1, param[3] + 1};
  // int nc[3] = {param[4], param[5], param[6]};

  // number of glops in subcell
  int npgc = nx[0] * nx[1] * nx[2]; 
  // number of glops in macrocell:
  int npg = nc[0] * nc[1] * nc[2] * npgc; 
  
  // index in subcell: 
  int ipgc = ix[0] + nx[0] * (ix[1] + nx[1] * ix[2]); 
  // index of subcell in macrocell:
  int nsubcell = ic[0] + nc[0] * (ic[1] + nc[1] * ic[2]);

  // index of point in macrocell:
  int ipg = ipgc + npgc * nsubcell; 
  return iv + m * (ipg + npg * elem);
}
#pragma end_opencl

// Given a the index ipg of a poing in a subcell, determine the three
// logical coordinates of that point in the subcell.
/* void ipg_to_xyz(int ipg, int *p, int *npg) */
/* { */
/*   p[0] = ipg % npg[0]; */
/*   p[1] = (ipg / npg[0]) % npg[1]; */
/*   p[2] = ipg / npg[0] / npg[1]; */
/* }  */

real min_grid_spacing(field *f)
{
  real hmin = FLT_MAX;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    real vol = 0, surf = 0;
    // Get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // Loop on the glops (for numerical integration)
    for(int ipg = 0; ipg < NPG(f->deg, f->raf);
	ipg++) {
      real xpgref[3], wpg;
      // Get the coordinates of the Gauss point
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      real codtau[3][3], dtau[3][3];
      Ref2Phy(physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      NULL, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      real det = dot_product(dtau[0], codtau[0]);
      vol += wpg * det;
    }
    for(int ifa = 0; ifa < 6; ifa++) {
      // loop on the faces
      for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, ifa); ipgf++) {
	real xpgref[3], wpg;
	// get the coordinates of the Gauss point
	ref_pg_face(f->deg, f->raf, ifa, ipgf, xpgref, &wpg, NULL);
	real vnds[3];
	{
	  real codtau[3][3], dtau[3][3];
	  Ref2Phy(physnode,
		  xpgref,
		  NULL, ifa, // dpsiref, ifa
		  NULL, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
	}
	surf += norm(vnds) * wpg;
      }
    }
    hmin = hmin < vol/surf ? hmin : vol/surf;

  }
 
  // Now take into account the polynomial degree and the refinement
  int maxd = f->interp_param[1];
  maxd = maxd > f->interp_param[2] ? maxd : f->interp_param[2];
  maxd = maxd > f->interp_param[3] ? maxd : f->interp_param[3];

  hmin /= ((maxd + 1) * f->interp_param[4]);

  return hmin;
}

void init_empty_field(field *f)
{
  f->period[0] = -1;
  f->period[1] = -1;
  f->period[2] = -1;
#ifdef _WITH_OPENCL
  f->use_source_cl = false;
#endif
}

void init_data(field *f)
{
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      real xpg[3];
      real xref[3], omega;
      ref_pg_vol(f->deg, f->raf, ipg, xref, &omega, NULL);
      real dtau[3][3];
      Ref2Phy(physnode,
	      xref,
	      0, -1, // dphiref, ifa
              xpg, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds
      { // Check the reverse transform at all the GLOPS
 	real xref2[3];
	Phy2Ref(physnode, xpg, xref2);
	assert(Dist(xref, xref2) < 1e-8);
      }

      real w[f->model.m];
      f->model.InitData(xpg, w);
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ie, ipg, iv);
	f->wn[imem] = w[iv];
      }
    }
  }
}

#ifdef _WITH_OPENCL
void set_physnodes_cl(field *f) 
{
  const int nmacro = f->macromesh.nbelems;
  real buf_size = sizeof(real) * 60 * nmacro;
  real *physnode = malloc(buf_size);
  
  for(int ie = 0; ie < nmacro; ++ie) {
    int ie20 = 20 * ie;
    for(int inoloc = 0; inoloc < 20; ++inoloc) {
      int ino = 3 * f->macromesh.elem2node[ie20 + inoloc];
      real *iphysnode = physnode + 3 * ie20 + 3 * inoloc;
      real *nodeino = f->macromesh.node + ino;
      iphysnode[0] = nodeino[0];
      iphysnode[1] = nodeino[1];
      iphysnode[2] = nodeino[2];
    }
  }
  
  cl_int status;
  status = clEnqueueWriteBuffer(f->cli.commandqueue,
  				f->physnodes_cl, // cl_mem buffer,
  				CL_TRUE,// cl_bool blocking_read,
  				0, // size_t offset
  				buf_size, // size_t cb
  				physnode, //  	void *ptr,
  				0, 0, 0);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  free(physnode);
}
#endif


void init_field_cl(field *f)
{
  printf("init_field_cl\n");
  InitCLInfo(&f->cli, nplatform_cl, ndevice_cl);
  cl_int status;

  f->wn_cl = clCreateBuffer(f->cli.context,
			    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
			    sizeof(real) * f->wsize,
			    f->wn,
			    &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dtwn_cl = clCreateBuffer(f->cli.context,
			      CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
			      sizeof(real) * f->wsize,
			      f->dtwn,
			      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->param_cl = clCreateBuffer(f->cli.context,
			       CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			       sizeof(int) * 7,
			       f->interp_param,
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Allocate and fill buffer for all macrocell geometries.
  const int nmacro = f->macromesh.nbelems;
  const size_t buf_size = sizeof(real) * 60 * nmacro;
  f->physnodes_cl = clCreateBuffer(f->cli.context,
				   CL_MEM_READ_ONLY,
				   buf_size,
				   NULL,
				   &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  set_physnodes_cl(f);

  // Program compilation
  char *strprog;
  GetOpenCLCode();
  ReadFile("schnaps.cl", &strprog);

  printf("\t%s\n", numflux_cl_name);
  //printf("\t%s\n", strprog);

  // If the source term is set (via set_source_CL) then add it to the
  // buildoptions and compile using the new buildoptions.
  if(f->use_source_cl) {
    char *temp;
    int len0 = strlen(cl_buildoptions);
    char *D_SOURCE_FUNC = " -D_SOURCE_FUNC=";
    int len1 = strlen(D_SOURCE_FUNC);
    int len2 = strlen(f->sourcename_cl);
    temp = calloc(sizeof(char), len0 + len1 + len2 + 2);
    printf("lens=%d %d %d \n",len0,len1,len2);
    strcat(temp, cl_buildoptions);
    strcat(temp, D_SOURCE_FUNC);
    strcat(temp, f->sourcename_cl);
    strcat(temp, " ");
    BuildKernels(&f->cli, strprog, temp);
  } else {
    printf("No source term\n");
    BuildKernels(&f->cli, strprog, cl_buildoptions);
  }

  f->dummykernel = clCreateKernel(f->cli.program,
				  "dummy",
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  
  f->dgmass = clCreateKernel(f->cli.program,
			     "DGMass",
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  f->dgflux = clCreateKernel(f->cli.program,
			     "DGFlux",
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dgvolume = clCreateKernel(f->cli.program,
			       "DGVolume",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dgsource = clCreateKernel(f->cli.program,
			       "DGSource",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dginterface = clCreateKernel(f->cli.program,
				  "DGMacroCellInterface",
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dgboundary = clCreateKernel(f->cli.program,
				 "DGBoundary",
				 &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->RK_out_CL = clCreateKernel(f->cli.program,
				"RK_out_CL",
				&status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->RK4_final_stage = clCreateKernel(f->cli.program,
				      "RK4_final_stage",
				      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->RK_in_CL = clCreateKernel(f->cli.program,
			       "RK_in_CL",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->zero_buf = clCreateKernel(f->cli.program,
			       "set_buffer_to_zero",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  
  // Initialize events. // FIXME: free on exit
  //f->clv_zbuf = clCreateUserEvent(f->cli.context, &status);
  
  const int ninterfaces = f->macromesh.nmacrointerfaces;
  if(ninterfaces > 0) {
    f->clv_mci = calloc(ninterfaces, sizeof(cl_event));
    for(int ifa = 0; ifa < ninterfaces; ++ifa) {
      //f->clv_mci[ifa] = clCreateUserEvent(f->cli.context, &status);
    }
  }
    
  const int nbound = f->macromesh.nboundaryfaces;
  if(nbound > 0) {
    f->clv_boundary = calloc(nbound, sizeof(cl_event));
    for(int ifa = 0; ifa < nbound; ++ifa) {
      //f->clv_boundary[ifa] = clCreateUserEvent(f->cli.context, &status);
    }
  }
  
  f->clv_mass = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    //f->clv_mass[ie] = clCreateUserEvent(f->cli.context, &status);
  }
    
  f->clv_flux0 = calloc(nmacro, sizeof(cl_event));
  f->clv_flux1 = calloc(nmacro, sizeof(cl_event));
  f->clv_flux2 = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    /* f->clv_flux0[ie] = clCreateUserEvent(f->cli.context, &status); */
    /* f->clv_flux1[ie] = clCreateUserEvent(f->cli.context, &status); */
    /* f->clv_flux2[ie] = clCreateUserEvent(f->cli.context, &status); */
  }
  
  f->clv_volume = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    //f->clv_volume[ie] = clCreateUserEvent(f->cli.context, &status);
  }

  f->clv_source = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    //f->clv_source[ie] = clCreateUserEvent(f->cli.context, &status);
  }
    
  // Set timers to zero
  f->zbuf_time = 0;
  f->mass_time = 0;
  f->vol_time = 0;
  f->flux_time = 0;
  f->minter_time = 0;
  f->boundary_time = 0;
  f->source_time = 0;
  f->rk_time = 0;

  // Set roofline counts to zero
  f->flops_vol = 0;
  f->flops_flux = 0;
  f->flops_mass = 0; 
  f->reads_vol = 0;
  f->reads_flux = 0;
  f->reads_mass = 0; 

  printf("... init_field_cl done\n");
}

void Initfield(field *f)
{
  int nmem = f->model.m * f->macromesh.nbelems * NPG(f->deg, f->raf);
  f->wsize = nmem;

  double g_memsize = nmem * sizeof(real) * 1e-9;
  if(sizeof(real) == sizeof(double))
    printf("Allocating %d doubles per array (%f GB).\n", nmem, g_memsize);
  else
    printf("Allocating %d floats per array (%f GB)\n", nmem, g_memsize);

  f->wn = calloc(nmem, sizeof(real));
  assert(f->wn);
  f->dtwn = calloc(nmem, sizeof(real));
  assert(f->dtwn);
  f->Diagnostics = NULL;
  f->pre_dtfield = NULL;
  f->post_dtfield = NULL;
  f->update_after_rk = NULL;
  f->model.Source = NULL;
  f->pic = NULL;

  // TODO: move this to the integrator code
  f->tnow=0;
  f->itermax=0;
  f->iter_time=0;
  f->nb_diags = 0;

  f->tnow = 0;

  init_data(f);

  // Compute cfl parameter min_i vol_i/surf_i
  f->hmin = min_grid_spacing(f);

  printf("hmin=%f\n", f->hmin);

  // Allocate and set MacroFaces
  f->mface = calloc(f->macromesh.nbfaces, sizeof(MacroFace*));
  for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++) {
    f->mface[ifa].first = ifa;
    f->mface[ifa].last_p1 = ifa + 1;
  }

  // Allocate and set MacroCells
  f->mcell = calloc(f->macromesh.nbelems, sizeof(MacroCell*));
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    f->mcell[ie].first = ie;
    f->mcell[ie].last_p1 = ie + 1;
  }

#ifdef _WITH_OPENCL
  // opencl inits
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable; OpenCL initialization disabled.\n");
  } else {

    init_field_cl(f);
  }
#endif // _WITH_OPENCL
  
  printf("field init done\n");
}

// This is the destructor for a field
void free_field(field *f) 
{
  free(f->mcell);
  free(f->mface);

#ifdef _WITH_OPENCL
  //cl_int status;

  /* status = clReleaseMemObject(f->physnode_cl); */
  /* if(status < CL_SUCCESS) printf("%s\n", clErrorString(status)); */
  /* assert(status >= CL_SUCCESS); */
#endif
}

// Display the field on screen
void Displayfield(field *f)
{
  printf("Display field...\n");
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    printf("elem %d\n", ie);
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      real xref[3], wpg;
      ref_pg_vol(f->deg, f->raf, ipg, xref, &wpg, NULL);

      printf("Gauss point %d %f %f %f \n", ipg, xref[0], xref[1], xref[2]);
      printf("dtw= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf,f->model.m,  ie, ipg, iv);
	printf("%f ", f->dtwn[imem]);
      }
      printf("\n");
      printf("w= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf,f->model.m,  ie, ipg, iv);
	printf("%f ", f->wn[imem]);
      }
      printf("\n");
    }
  }
}

// Save the results in a text file
// in order plot it with Gnuplot
void Gnuplot(field* f,int dir, real fixval, char* filename)
{
  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  printf("Save for Gnuplot...\n");
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {

    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {

      real xref[3], xphy[3], wpg;
      real dtau[3][3];
      ref_pg_vol(f->deg, f->raf, ipg, xref, &wpg, NULL);

      Ref2Phy(physnode,
	      xref,
	      0, -1, // dphiref, ifa
	      xphy, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds

      if(xphy[dir] > -(fixval + 0.0001) && xphy[dir] < (fixval + 0.00001)){

	fprintf(gmshfile, "%f ",xphy[1-dir]);

	for(int iv = 0; iv < f->model.m; iv++) {
	  int imem = f->varindex(f->deg, f->raf,f->model.m,  ie, ipg, iv);
	  fprintf(gmshfile, "%f ",f->wn[imem]);
	}
	fprintf(gmshfile, "\n");

      }
    }
  }
  fclose(gmshfile);
}

// Save the results in the gmsh format typplot: index of the plotted
// variable int compare == true -> compare with the exact value.  If
// fieldname is NULL, then the fieldname is typpplot.
void Plotfield(int typplot, int compare, field* f, char *fieldname,
	       char *filename) {

  real hexa64ref[3 * 64] = { 
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;

  int *elem2nodes = f->macromesh.elem2node;
  real *node = f->macromesh.node;

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int raf[3] = {f->interp_param[4], f->interp_param[5], f->interp_param[6]};
  // Refinement size in each direction
  real hh[3] = {1.0 / raf[0], 1.0 / raf[1], 1.0 / raf[2]};

  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(real));

  int nb_plotnodes = f->macromesh.nbelems * raf[0] * raf[1] * raf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);

  real *value = malloc(nb_plotnodes * sizeof(real));
  assert(value);
  int nodecount = 0;

  // Nodes
  int npgv = NPG(f->deg, f->raf);
  for(int i = 0; i < f->macromesh.nbelems; i++) {
    // Get the nodes of element L
    int nnodes = 20;
    real physnode[nnodes][3];
    for(int ino = 0; ino < nnodes; ino++) {
      int numnoe = elem2nodes[nnodes * i + ino];
      for(int ii = 0; ii < 3; ii++) {
        physnode[ino][ii] = node[3 * numnoe + ii];
      }
    }

    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < raf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < raf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < raf[2]; icL[2]++) {

	  for(int ino = 0; ino < 64; ino++) {
	    real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
			   hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
			   hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
	    
	    for(int ii = 0; ii < 3; ii++) {
	      assert(Xr[ii] < 1 +  1e-10);
	      assert(Xr[ii] > -1e-10);
	    }

	    real Xphy[3];
	    Ref2Phy(physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

	    real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};

	    value[nodecount] = 0;
	    real testpsi = 0;
	    for(int ib = 0; ib < npgv; ib++) {
	      real psi;
	      psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
	      testpsi += psi;
	      int vi = f->varindex(f->deg, f->raf, f->model.m ,i, ib, typplot);
	      value[nodecount] += psi * f->wn[vi];
	    }
	    assert(fabs(testpsi-1) < 1e-10);

	    // Compare with an exact solution
	    if (compare) {
	      real wex[f->model.m];
	      f->model.ImposedData(Xphy, f->tnow, wex);
	      value[nodecount] -= wex[typplot];

	    }
	    nodecount++;
	    fprintf(gmshfile, "%d %f %f %f\n", nodecount,
		    Xplot[0], Xplot[1], Xplot[2]);
	  }
	}
      }
    }
  }

  fprintf(gmshfile, "$EndNodes\n");

  // Elements
  fprintf(gmshfile, "$Elements\n");
  fprintf(gmshfile, "%d\n", 
	  f->macromesh.nbelems * raf[0] * raf[1] * raf[2]);

  int elm_type = 92;
  int num_tags = 0;

  // fwrite((char*) &elm_type, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_elm_follow, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_tags, sizeof(int), 1, gmshfile);

  for(int i = 0; i < f->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < raf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < raf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < raf[2]; icL[2]++) {
	  // Get the subcell id
	  int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);

	  // Global subcell id
	  int numelem = ncL + i * raf[0] * raf[1] * raf[2] + 1;

	  //fwrite((char*) &numelem, sizeof(int), 1, gmshfile);
	  fprintf(gmshfile, "%d ", numelem);
	  fprintf(gmshfile, "%d ", elm_type);
	  fprintf(gmshfile, "%d ", num_tags);

	  for(int ii = 0; ii < 64; ii++) {
	    int numnoe = 64 * (i * raf[0] * raf[1] * raf[2] + ncL) + ii  + 1;
	    //fwrite((char*) &numnoe, sizeof(int), 1, gmshfile);
	    fprintf(gmshfile, "%d ", numnoe);
	  }
	  fprintf(gmshfile, "\n");
	}
      }
    }
  }

  fprintf(gmshfile, "$EndElements\n");

  // Now display data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", typplot);
  else 
    fprintf(gmshfile, "\"field: %s\"\n", fieldname);

  real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);

  for(int ino = 0; ino < nb_plotnodes; ino++) {
    //fwrite(const void *ptr, size_t size_of_elements,
    // size_t number_of_elements, FILE *a_file);
    //fwrite((char*) &nodenumber, sizeof(int), 1, gmshfile);
    //fwrite((char*) &value, sizeof(real), 1, gmshfile);
    //fprintf(gmshfile, "%d %f\n", nodenumber, value);
    fprintf(gmshfile, "%d %f\n", ino + 1, value[ino]);
  }

  fprintf(gmshfile, "\n$EndNodeData\n");

  fclose(gmshfile);
  free(value);
}

// Compute inter-subcell fluxes
void DGSubCellInterface(int ie, field *f, real *w, real *dtw) 
{
  // get the physical nodes of element ie
  real physnode[20][3];
  for(int inoloc = 0; inoloc < 20; inoloc++) {
    int ino = f->macromesh.elem2node[20 * ie + inoloc];
    physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
    physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
    physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
  }

  const int raf[3] = {f->interp_param[4],
		      f->interp_param[5],
		      f->interp_param[6]};
  const int deg[3] = {f->interp_param[1],
		      f->interp_param[2],
		      f->interp_param[3]};
  const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
  const int m = f->model.m;

  // Loop on the subcells
  for(int icL0 = 0; icL0 < raf[0]; icL0++) {
    for(int icL1 = 0; icL1 < raf[1]; icL1++) {
      for(int icL2 = 0; icL2 < raf[2]; icL2++) {

	int icL[3] = {icL0, icL1, icL2};

	// Get the left subcell id
	int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);
	// First glop index in the subcell
	int offsetL = npg[0] * npg[1] * npg[2] * ncL;

	// Sweeping subcell faces in the three directions
	for(int dim0 = 0; dim0 < 3; dim0++) { 
	    
	  // Compute the subface flux only if we do not touch the
	  // subcell boundary along the current direction dim0
	  if (icL[dim0] != raf[dim0] - 1) {
	    int icR[3] = {icL[0], icL[1], icL[2]};
	    // The right cell index corresponds to an increment in
	    // the dim0 direction
	    icR[dim0]++;
	    int ncR = icR[0] + raf[0] * (icR[1] + raf[1] * icR[2]);
	    int offsetR = npg[0] * npg[1] * npg[2] * ncR;

	    // FIXME: write only write to L-values (and do both
	    // faces) to parallelise better.

	    const int altdim1[3] = {1, 0, 0};
	    const int altdim2[3] = {2, 2, 1};

	    // now loop on the left glops of the subface
	    //int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;
	    int dim1 = altdim1[dim0];
	    int dim2 = altdim2[dim0];
	    int iL[3];
	    iL[dim0] = deg[dim0];
	    for(iL[dim2] = 0; iL[dim2] < npg[dim2]; iL[dim2]++) {
	      for(iL[dim1] = 0; iL[dim1] < npg[dim1]; iL[dim1]++) {
		// find the right and left glops volume indices

		int iR[3] = {iL[0], iL[1], iL[2]};
		iR[dim0] = 0;

		int ipgL = offsetL 
		  + iL[0] + (deg[0] + 1) * (iL[1] + (deg[1] + 1) * iL[2]);
		int ipgR = offsetR 
		  + iR[0] + (deg[0] + 1) * (iR[1] + (deg[1] + 1) * iR[2]);
		//printf("ipgL=%d ipgR=%d\n", ipgL, ipgR);

		// Compute the normal vector for integrating on the
		// face
		real vnds[3];
		{
		  real xref[3], wpg3;
		  ref_pg_vol(f->deg, f->raf, ipgL, xref, &wpg3, NULL);
		  // mapping from the ref glop to the physical glop
		  real dtau[3][3], codtau[3][3];
		  Ref2Phy(physnode,
			  xref,
			  NULL, // dphiref
			  -1,  // ifa
			  NULL, // xphy
			  dtau,
			  codtau,
			  NULL, // dphi
			  NULL);  // vnds
		  // we compute ourself the normal vector because we
		  // have to take into account the subcell surface

		  real h1h2 = 1. / raf[dim1] / raf[dim2];
		  vnds[0] = codtau[0][dim0] * h1h2;
		  vnds[1] = codtau[1][dim0] * h1h2;
		  vnds[2] = codtau[2][dim0] * h1h2;
		}

		// numerical flux from the left and right state and
		// normal vector
		real wL[m], wR[m], flux[m];
		for(int iv = 0; iv < m; iv++) {
		  // TODO change the varindex signature
		  int imemL = f->varindex(f->deg, f->raf, f->model.m,
					  ie, ipgL, iv); 
		  int imemR = f->varindex(f->deg, f->raf, f->model.m,
					  ie, ipgR, iv);
		  // end TODO
		  wL[iv] = w[imemL];
		  wR[iv] = w[imemR];
		}
		f->model.NumFlux(wL, wR, vnds, flux);

		// subcell ref surface glop weight
		real wpg
		  = wglop(deg[dim1], iL[dim1])
		  * wglop(deg[dim2], iL[dim2]);

		/* printf("vnds %f %f %f flux %f wpg %f\n", */
		/* 	 vnds[0], vnds[1], vnds[2], */
		/* 	 flux[0], wpg); */

		// finally distribute the flux on the two sides
		for(int iv = 0; iv < m; iv++) {
		  // TO DO change the varindex signature
		  int imemL = f->varindex(f->deg, f->raf, f->model.m,
					  ie, ipgL, iv);
		  int imemR = f->varindex(f->deg, f->raf, f->model.m,
					  ie, ipgR, iv);
		  // end TO DO
		  dtw[imemL] -= flux[iv] * wpg;
		  dtw[imemR] += flux[iv] * wpg;
		}

	      }  // face yhat loop
	    } // face xhat loop
	  } // endif internal face
	} // dim loop
      } // subcell icl2 loop
    } // subcell icl1 loop
  } // subcell icl0 loop
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
void DGMacroCellInterface(void *mc, field *f, real *w, real *dtw, real tnow) 
{
  MacroFace *mface = (MacroFace*) mc;
  MacroMesh *msh = &f->macromesh;
  const unsigned int m = f->model.m;

  // Assembly of the surface terms loop on the macrocells faces
  for (int ifa = mface->first; ifa < mface->last_p1; ifa++) {
    int ieL = msh->face2elem[4 * ifa + 0];
    int locfaL = msh->face2elem[4 * ifa + 1];

    // Get the physical nodes of element ieL
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = msh->elem2node[20 * ieL + inoloc];
      physnode[inoloc][0] = msh->node[3 * ino + 0];
      physnode[inoloc][1] = msh->node[3 * ino + 1];
      physnode[inoloc][2] = msh->node[3 * ino + 2];
    }

    int ieR = msh->face2elem[4 * ifa + 2];
    int locfaR = msh->face2elem[4 * ifa + 3];
    real physnodeR[20][3];
    if (ieR >= 0) {
      for(int inoloc = 0; inoloc < 20; inoloc++) {
        int ino = msh->elem2node[20 * ieR + inoloc];
        physnodeR[inoloc][0] = msh->node[3 * ino + 0];
        physnodeR[inoloc][1] = msh->node[3 * ino + 1];
        physnodeR[inoloc][2] = msh->node[3 * ino + 2];
      }
    }

    // Loop over the points on a single macro cell interface.
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ipgfL = 0; ipgfL < NPGF(f->deg, f->raf, locfaL); ipgfL++) {
      real xpgref[3], xpgref_in[3], wpg;
      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      //ref_pg_face(f->deg, f->raf, locfaL, ipgfL, xpgref, &wpg, xpgref_in);

      int ipgL = ref_pg_face(f->deg, f->raf, locfaL, ipgfL, xpgref, &wpg,
			     xpgref_in);
      
      // Recover the volume gauss point from the face index
      real flux[m];
      real wL[m];

      // Normal vector at gauss point ipgL
      real vnds[3], xpg[3];
      {
	real dtau[3][3], codtau[3][3];
	Ref2Phy(physnode,
		xpgref,
		NULL, locfaL, // dpsiref, ifa
		xpg, dtau,
		codtau, NULL, vnds); // codtau, dpsi, vnds
      }

      if (ieR >= 0) {  // the right element exists
        real xrefL[3];
	{
	  real xpg_in[3];
	  Ref2Phy(physnode,
		  xpgref_in,
		  NULL, -1, // dpsiref, ifa
		  xpg_in, NULL,
		  NULL, NULL, NULL); // codtau, dpsi, vnds
	  PeriodicCorrection(xpg_in,f->macromesh.period);
	  Phy2Ref(physnodeR, xpg_in, xrefL);

	}

        int ipgR = ref_ipg(f->deg, f->raf, xrefL);

	real wR[m];
        for(int iv = 0; iv < m; iv++) {
	  int imemL = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv);
	  wL[iv] = w[imemL];
          int imemR = f->varindex(f->deg, f->raf, f->model.m, ieR, ipgR, iv);
          wR[iv] = w[imemR];
        }

        f->model.NumFlux(wL, wR, vnds, flux);

	// Add flux to both sides
	for(int iv = 0; iv < m; iv++) {
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv);
          int imemR = f->varindex(f->deg, f->raf, f->model.m, ieR, ipgR, iv);
	  dtw[imemL] -= flux[iv] * wpg;
          dtw[imemR] += flux[iv] * wpg;
	}

      } else { 
	// The point is on the boundary.
	for(int iv = 0; iv < m; iv++) {
	  int imemL = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv);
	  wL[iv] = w[imemL];
	}

        f->model.BoundaryFlux(xpg, tnow, wL, vnds, flux);

	for(int iv = 0; iv < m; iv++) {
	  int imemL = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv);
	  dtw[imemL] -= flux[iv] * wpg;
	}
      }

    }

  }
}

// Apply division by the mass matrix
void DGMass(int ie, field *f, real *dtw) 
{
  int m = f->model.m;

  // get the physical nodes of element ie
  real physnode[20][3];
  for(int inoloc = 0; inoloc < 20; inoloc++) {
    int ino = f->macromesh.elem2node[20 * ie + inoloc];
    physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
    physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
    physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
  }
  for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
    real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
    Ref2Phy(physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    real det = dot_product(dtau[0], codtau[0]);
    for(int iv = 0; iv < f->model.m; iv++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m, ie, ipg, iv);
      dtw[imem] /= (wpg * det);
    }
  }
}

// Apply the source term
void DGSource(int ie, field *f, real *w, real *dtw) 
{
  if (f->model.Source == NULL) {
    return;
  }
  
  const int m = f->model.m;

  // Get the physical nodes of element ie
  real physnode[20][3];
  for(int inoloc = 0; inoloc < 20; inoloc++) {
    int ino = f->macromesh.elem2node[20 * ie + inoloc];
    physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
    physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
    physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
  }
  for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
    real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
    Ref2Phy(physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    real wL[m], source[m];
    for(int iv = 0; iv < m; ++iv){
      int imem = f->varindex(f->deg, f->raf, f->model.m, ie, ipg, iv);
      wL[iv] = w[imem];
    }
      
    f->model.Source(xphy, f->tnow, wL, source);
      
    for(int iv = 0; iv < m; ++iv) {
      int imem = f->varindex(f->deg, f->raf, f->model.m, ie, ipg, iv);
      dtw[imem] += source[iv];
	
    }
  }

}

// Compute the Discontinuous Galerkin volume terms, fast version
void DGVolume(int ie, field *f, real *w, real *dtw) 
{
  // get the physical nodes of element ie
  real physnode[20][3];
  for(int inoloc = 0; inoloc < 20; inoloc++) {
    int ino = f->macromesh.elem2node[20*ie+inoloc];
    physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
    physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
    physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
  }

  const int m = f->model.m;
  const int deg[3] = {f->deg[0],
		      f->deg[1],
		      f->deg[2]};
  const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
  const int raf[3] = {f->raf[0],
		      f->raf[1],
		      f->raf[2]};

  const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Loop on the subcells
  for(int icL0 = 0; icL0 < raf[0]; icL0++) {
    for(int icL1 = 0; icL1 < raf[1]; icL1++) {
      for(int icL2 = 0; icL2 < raf[2]; icL2++) {

	int icL[3] = {icL0, icL1, icL2};
	// get the L subcell id
	int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);
	// first glop index in the subcell

	int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	
	// compute all of the xref for the subcell
	real *xref0 = malloc(sc_npg * sizeof(real));
	real *xref1 = malloc(sc_npg * sizeof(real));
	real *xref2 = malloc(sc_npg * sizeof(real));
	real *omega = malloc(sc_npg * sizeof(real));
	int *imems = malloc(m * sc_npg * sizeof(int));
	int pos = 0;
	for(unsigned int p = 0; p < sc_npg; ++p) {
	  real xref[3];
	  real tomega;

	  ref_pg_vol(f->deg, f->raf, offsetL + p, xref, &tomega, NULL);
	  xref0[p] = xref[0];
	  xref1[p] = xref[1];
	  xref2[p] = xref[2];
	  omega[p] = tomega;

	  for(int im = 0; im < m; ++im) {
	    imems[pos++] = f->varindex(f->deg, f->raf, f->model.m,
				       ie, offsetL + p, im);
	  }
	}

	// loop in the "cross" in the three directions
	for(int dim0 = 0; dim0 < 3; dim0++) {
	  // point p at which we compute the flux

	  for(int p0 = 0; p0 < npg[0]; p0++) {
	    for(int p1 = 0; p1 < npg[1]; p1++) {
	      for(int p2 = 0; p2 < npg[2]; p2++) {
		real wL[m], flux[m];
		int p[3] = {p0, p1, p2};
		int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
		for(int iv = 0; iv < m; iv++) {
		  int imemL = f->varindex(f->deg, f->raf, f->model.m, 
					  ie, ipgL, iv);
		  //wL[iv] = w[imems[m * (ipgL - offsetL) + iv]];
		}
		int q[3] = {p[0], p[1], p[2]};
		// loop on the direction dim0 on the "cross"
		for(int iq = 0; iq < npg[dim0]; iq++) {
		  q[dim0] = (p[dim0] + iq) % npg[dim0];
		  real dphiref[3] = {0, 0, 0};
		  // compute grad phi_q at glop p
		  dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) 
		    * raf[dim0];

		  real xrefL[3] = {xref0[ipgL - offsetL],
				   xref1[ipgL - offsetL],
				   xref2[ipgL - offsetL]};
		  real wpgL = omega[ipgL - offsetL];
		  /* real xrefL[3], wpgL; */
		  /* ref_pg_vol(f->interp_param+1,ipgL,xrefL, &wpgL, NULL); */

		  // mapping from the ref glop to the physical glop
		  real dtau[3][3], codtau[3][3], dphiL[3];
		  Ref2Phy(physnode,
			  xrefL,
			  dphiref, // dphiref
			  -1,  // ifa
			  NULL, // xphy
			  dtau,
			  codtau,
			  dphiL, // dphi
			  NULL);  // vnds

		  f->model.NumFlux(wL, wL, dphiL, flux);

		  int ipgR = offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		  for(int iv = 0; iv < m; iv++) {
		    int imemR = f->varindex(f->deg, f->raf, f->model.m,
					    ie, ipgR, iv);
		    int temp = m * (ipgR - offsetL) + iv;  
		    assert(imemR == imems[temp]);
		    dtw[imems[temp]] += flux[iv] * wpgL;
		  }
		} // iq
	      } // p2
	    } // p1
	  } // p0

	} // dim loop

	free(omega);
	free(xref0);
	free(xref1);
	free(xref2);
	free(imems);

      } // icl2
    } //icl1
  } // icl0
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field
void dtfield(field *f, real *w, real *dtw, real tnow)
{
  if(f->pre_dtfield != NULL) {
    f->pre_dtfield(f, w);
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for(int iw = 0; iw < f->wsize; iw++) {
    dtw[iw] = 0;
  }

  for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++) {
    DGMacroCellInterface((void*) (f->mface + ifa), f, w, dtw, tnow);
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    DGSubCellInterface(ie, f, w, dtw);
    DGVolume(ie, f, w, dtw);
    DGMass(ie, f, dtw);
    DGSource(ie, f, w, dtw);
  }

  if(f->post_dtfield != NULL) {
    f->post_dtfield(f, w);
  }
}

// An out-of-place RK step
void RK_out(real *dest, real *fwn, real *fdtwn, const real dt, 
	    const int sizew)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < sizew; iw++) {
    dest[iw] = fwn[iw] + dt * fdtwn[iw];
  }
}

// An in-place RK step
void RK_in(real *fwnp1, real *fdtwn, const real dt, const int sizew)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < sizew; iw++) {
    fwnp1[iw] += dt * fdtwn[iw];
  }
}

real set_dt(field *f)
{
  return f->model.cfl * f->hmin / f->vmax; 
}

// Time integration by a second-order Runge-Kutta algorithm
void RK2(field *f, real tmax, real dt) 
{
  if(dt <= 0)
    dt = set_dt(f);

  real tnow = 0.0;

  f->itermax = tmax / dt;
  int size_diags;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(f->deg, f->raf);
  int iter = 0;

  real *wnp1 = calloc(f->wsize, sizeof(real));
  assert(wnp1);

  // FIXME: remove
  size_diags = f->nb_diags * f->itermax;
  f->iter_time = iter;

  if(f->nb_diags != 0)
    f->Diagnostics = malloc(size_diags * sizeof(real));

  while(tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", tnow, iter, f->itermax, dt);

    dtfield(f, f->wn, f->dtwn, tnow);
    RK_out(wnp1, f->wn, f->dtwn, 0.5 * dt, sizew);

    tnow += 0.5 * dt;

    dtfield(f, wnp1, f->dtwn, tnow);
    RK_in(f->wn, f->dtwn, dt, sizew);

    tnow += 0.5 * dt;

    if(f->update_after_rk != NULL)
      f->update_after_rk(f, f->wn);

    iter++;
    f->iter_time=iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", tnow, iter, f->itermax, dt);
  free(wnp1);
}

void RK4_final_inplace(real *w, real *l1, real *l2, real *l3, 
		       real *dtw, const real dt, const int sizew)
{
  const real b = -1.0 / 3.0;
  const real a[] = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, dt / 6.0};
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < sizew; ++i) {
    w[i] = 
      b * w[i] +
      a[0] * l1[i] +
      a[1] * l2[i] +
      a[2] * l3[i] +
      a[3] * dtw[i];
  }
}

// Time integration by a fourth-order Runge-Kutta algorithm
void RK4(field *f, real tmax, real dt) 
{
  if(dt <= 0)
    dt = set_dt(f);

  real tnow = 0.0;

  f->itermax = tmax / dt;
  int size_diags;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(f->deg, f->raf);
  int iter = 0;

  // Allocate memory for RK time-stepping
  real *l1, *l2, *l3;
  l1 = calloc(sizew, sizeof(real));
  l2 = calloc(sizew, sizeof(real));
  l3 = calloc(sizew, sizeof(real));
  
  size_diags = f->nb_diags * f->itermax;
  f->iter_time = iter;
  
  if(f->nb_diags != 0)
    f->Diagnostics = malloc(size_diags * sizeof(real));
  
  while(tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", tnow, iter, f->itermax, dt);

    // l_1 = w_n + 0.5dt * S(w_n, t_0)
    dtfield(f, f->wn, f->dtwn, tnow);
    RK_out(l1, f->wn, f->dtwn, 0.5 * dt, sizew);

    tnow += 0.5 * dt;

    // l_2 = w_n + 0.5dt * S(l_1, t_0 + 0.5 * dt)
    dtfield(f, l1, f->dtwn, tnow);
    RK_out(l2, f->wn, f->dtwn, 0.5 * dt, sizew);

    // l_3 = w_n + dt * S(l_2, t_0 + 0.5 * dt)
    dtfield(f, l2, f->dtwn, tnow);
    RK_out(l3, f->wn, f->dtwn, dt, sizew);

    tnow += 0.5 * dt;

    // Compute S(l_3, t_0 + dt)
    dtfield(f, l3, f->dtwn, tnow);
    RK4_final_inplace(f->wn, l1, l2, l3, f->dtwn, dt, sizew);

    
    if(f->update_after_rk != NULL)
      f->update_after_rk(f, f->wn);
    
    iter++;
    f->iter_time=iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", tnow, iter, f->itermax, dt);

  free(l3);
  free(l2);
  free(l1);
}

// Compute the normalized L2 distance with the imposed data
real L2error(field *f) {
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  real error = 0;
  real mean = 0;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    // Get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ie, ipg, iv);
	w[iv] = f->wn[imem];
      }

      real wex[f->model.m];
      real wpg, det;
      // Compute wpg, det, and the exact solution
      { 
	real xphy[3], xpgref[3];
	real dtau[3][3], codtau[3][3];
	// Get the coordinates of the Gauss point
	ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
	Ref2Phy(physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	det = dot_product(dtau[0], codtau[0]);

	// Get the exact value
	f->model.ImposedData(xphy, f->tnow, wex);
      }

      for(int iv = 0; iv < f->model.m; iv++) {
	real diff = w[iv] - wex[iv];
	error += diff * diff * wpg * det;
	mean += wex[iv] * wex[iv] * wpg * det;
      }
    }
  }
  return sqrt(error) / (sqrt(mean)  + 1e-16);
}


// Compute the normalized L2 distance with the imposed data
real L2error_onefield(field *f, int nbfield) 
{
  real error = 0;
  real mean = 0;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    // Get the physical nodes of element ie

    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ie, ipg, iv);
	w[iv] = f->wn[imem];
      }

      real wex[f->model.m];
      real wpg, det;
      // Compute wpg, det, and the exact solution
      { 
	real xphy[3], xpgref[3];
	real dtau[3][3], codtau[3][3];
	// Get the coordinates of the Gauss point
	ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);

	real physnode[20][3];
	for(int inoloc = 0; inoloc < 20; inoloc++) {
	  int ino = f->macromesh.elem2node[20*ie+inoloc];
	  physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
	  physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
	  physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
	}

	Ref2Phy(physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	det = dot_product(dtau[0], codtau[0]);

	// Get the exact value
	f->model.ImposedData(xphy, f->tnow, wex);
      }

      int iv = nbfield;
      real diff = w[iv] - wex[iv];
      error += diff * diff * wpg * det;
      mean += w[iv] * w[iv] * wpg * det;
    }
  }
  return sqrt(error);
}


void InterpField(field *f, int ie, real *xref, real *w){
  const int raf[3] = {f->interp_param[4],
		      f->interp_param[5],
		      f->interp_param[6]};
  const int deg[3] = {f->interp_param[1],
		      f->interp_param[2],
		      f->interp_param[3]};

  for(int iv = 0; iv < f->model.m; iv++)
    w[iv] = 0;

  int is[3];

  for(int ii = 0; ii < 3; ii++){
    is[ii] = xref[ii] * raf[ii];
    assert(is[ii] < raf[ii] && is[ii]>= 0);
  }
  
  int npgv = NPG(f->deg, f->raf);
  // TODO: loop only on non zero basis function
  for(int ib = 0; ib < npgv; ib++) { 
    real psi;
    psi_ref_subcell(f->deg, f->raf, is, ib, xref, &psi, NULL);
    
    for(int iv=0;iv<f->model.m;iv++){
      int imem = f->varindex(f->deg, f->raf, f->model.m, ie, ib, iv);
      w[iv] += psi * f->wn[imem];
    }
  }
}
