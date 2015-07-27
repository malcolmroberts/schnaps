#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>




void InitSimulation(Simulation *simu, MacroMesh *mesh,
		    int *deg, int *raf, Model *model){

  simu->macromesh = *mesh;

  simu->tnow = 0;

  simu->fd = malloc(mesh->nbelems * sizeof(field));

  field *fd = simu->fd;

  int field_size = NPG(deg,raf) * model->m;

  simu->wsize = field_size * mesh->nbelems;

  printf("field_size = %d simusize = %d\n",field_size,simu->wsize);

  simu->w = calloc(simu->wsize, sizeof(real));
  simu->dtw = calloc(simu->wsize, sizeof(real));

  real *w = simu->w;
  real *dtw = simu->dtw;
  
  real physnode[20][3];
  
  simu->hmin = FLT_MAX;

  for(int ie=0; ie < mesh->nbelems; ++ie){
    for(int inoloc = 0; inoloc < 20; ++inoloc){
      int ino = mesh->elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = mesh->node[3 * ino + 0];
      physnode[inoloc][1] = mesh->node[3 * ino + 1];
      physnode[inoloc][2] = mesh->node[3 * ino + 2];
    }

    init_empty_field(fd + ie);

    fd[ie].period[0] = simu->macromesh.period[0];
    fd[ie].period[1] = simu->macromesh.period[1];
    fd[ie].period[2] = simu->macromesh.period[2];

    Initfield(fd + ie, *model, physnode, deg, raf,
	      w + ie * field_size, dtw + ie * field_size);

    simu->hmin = simu->hmin > fd[ie].hmin ? fd[ie].hmin : simu->hmin;
  }

  simu->pre_dtfields = NULL;

  simu->nb_diags = 0;

}


void DisplaySimulation(Simulation *simu){

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie){
    printf("Field %d\n",ie);

    Displayfield(simu->fd + ie);
  }

}

// Save the results in the gmsh format typplot: index of the plotted
// variable int compare == true -> compare with the exact value.  If
// fieldname is NULL, then the fieldname is typpplot.
void PlotFields(int typplot, int compare, Simulation* simu, char *fieldname,
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

  int *elem2nodes = simu->macromesh.elem2node;
  real *node = simu->macromesh.node;

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1], 
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1], 
		simu->fd[0].deg[2]};
  // Refinement size in each direction
  real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};

  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(real));

  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);

  real *value = malloc(nb_plotnodes * sizeof(real));
  assert(value);
  int nodecount = 0;

  // Nodes
  int npgv = NPG(deg, nraf);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    int nnodes = 20;

    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {

	  for(int ino = 0; ino < 64; ino++) {
	    real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
			     hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
			     hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
	    
	    for(int ii = 0; ii < 3; ii++) {
	      assert(Xr[ii] < 1 +  1e-10);
	      assert(Xr[ii] > -1e-10);
	    }

	    real Xphy[3];
	    Ref2Phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

	    real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};

	    value[nodecount] = 0;
	    real testpsi = 0;
	    for(int ib = 0; ib < npgv; ib++) {
	      real psi;
	      psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
	      testpsi += psi;
	      int vi = f->varindex(f->deg, f->raf, f->model.m, ib, typplot);
	      value[nodecount] += psi * f->wn[vi];
	    }
	    assert(fabs(testpsi-1) < _SMALL);

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
	  simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]);

  int elm_type = 92;
  int num_tags = 0;

  // fwrite((char*) &elm_type, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_elm_follow, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_tags, sizeof(int), 1, gmshfile);

  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
	  // Get the subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);

	  // Global subcell id
	  int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;

	  //fwrite((char*) &numelem, sizeof(int), 1, gmshfile);
	  fprintf(gmshfile, "%d ", numelem);
	  fprintf(gmshfile, "%d ", elm_type);
	  fprintf(gmshfile, "%d ", num_tags);

	  for(int ii = 0; ii < 64; ii++) {
	    int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
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

  /* for(int i = 0;i < f->macromesh.nbelems;i++) { */
  /*   for(int ino = 0;ino < 20;ino++) { */
  /* 	int numnoe = elem2nodes[nnodes*i+ino]; */
  /* 	for(int ii = 0;ii < 3;ii++) { */
  /* 	  physnode[ino][ii] = node[3 * numnoe+ii]; */
  /* 	} */
  /*   } */

  /*   // data at the eight nodes */
  /*   for(int ii = 0;ii < 64;ii++) { */
  /* 	int nodenumber = 64*i + ii  + 1; */

  /* 	Xr[0] = (real) (hexa64ref[3 * ii+0]) / 3; */
  /* 	Xr[1] = (real) (hexa64ref[3 * ii + 1]) / 3; */
  /* 	Xr[2] = (real) (hexa64ref[3 * ii+2]) / 3; */

  /* 	Ref2Phy(physnode, */
  /* 		Xr, */
  /* 		NULL, */
  /* 		-1, */
  /* 		Xphy, */
  /* 		NULL, */
  /* 		NULL, */
  /* 		NULL, */
  /* 		NULL); */


  /* 	real value = 0; */
  /* 	for(int ib = 0;ib < npgv;ib++) { */
  /* 	  real psi; */
  /* 	  psi_ref(f->deg, f->raf, ib, Xr, &psi, NULL); */

  /* 	  int vi  =  f->varindex(f->interp_param, i, ib, typplot); */
  /* 	  value += psi * f->wn[vi]; */
  /* 	} */

  /* 	// compare with an */
  /* 	// exact solution */
  /*     if (compare) { */
  /*       real wex[f->model.m]; */
  /*       f->model.ImposedData(Xphy, f->tnow, wex); */
  /*       value -= wex[typplot]; */
  /*     } */


  /* 	//fwrite(const void *ptr, size_t size_of_elements, */
  /* 	// size_t number_of_elements, FILE *a_file); */
  /* 	//fwrite((char*) &nodenumber, sizeof(int), 1, gmshfile); */
  /* 	//fwrite((char*) &value, sizeof(real), 1, gmshfile); */
  /* 	//fprintf(gmshfile, "%d %f\n", nodenumber, value); */
  /* 	fprintf(gmshfile, "%d %f\n", nodenumber, value); */
  /*   } */

  /* } */

  fprintf(gmshfile, "\n$EndNodeData\n");

  fclose(gmshfile);
  free(value);
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field
void DtFields(Simulation *simu, real *w, real *dtw) {

  if(simu->pre_dtfields != NULL) {
    simu->pre_dtfields(simu, w);
    //assert(1==2);
  }


#ifdef _OPENMP
#pragma omp parallel
#endif

  //real *w = simu->fd[0].wn;
  //real *dtw = simu->fd[0].dtwn;

  int fsize =  simu->wsize / simu->macromesh.nbelems;
  
  for(int iw = 0; iw < simu->wsize; iw++)
    dtw[iw] = 0;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  } 

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      fR = simu->fd + ieR;
      offsetR = fsize * ieR;
    }
         
    
 
    DGMacroCellInterface(locfaL,
			 fL, offsetL, fR, offsetR,
			 w, dtw);
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    DGSubCellInterface(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGVolume(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGMass(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGSource(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);

  }

  //if(f->post_dtfield != NULL) // FIXME: rename to after dtfield
    //f->post_dtfield(f, w);
}

real L2error(Simulation *simu) {

  real error = 0;
  real mean = 0;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {

    field *f = simu->fd + ie;

    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
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
	Ref2Phy(f->physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	det = dot_product(dtau[0], codtau[0]);

	// Get the exact value
	f->model.ImposedData(xphy, simu->tnow, wex);
      }

      for(int iv = 0; iv < f->model.m; iv++) {
	//for(int iv = 0; iv < 4; iv++) {   ///////error here for coil2d
	real diff = w[iv] - wex[iv];
       error += diff * diff * wpg * det;
        mean += w[iv] * w[iv] * wpg * det;
	//printf("ie=%d ipg=%d iv=%d err=%f \n",ie,ipg,iv,diff);
        }
    }
  }
  //printf("errl2=%f\n",sqrt(error) / (sqrt(mean)  + 1e-14));
  return sqrt(error) / (sqrt(mean)  + 1e-14);
}

// Time integration by a second-order Runge-Kutta algorithm
void RK2(Simulation *simu, real tmax){

  simu->dt = Get_Dt_RK(simu);


  real dt = simu->dt;

  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  real *wnp1 = calloc(simu->wsize, sizeof(real));
  assert(wnp1);

  // FIXME: remove
  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

  /* if(simu->nb_diags != 0) */
  /*   simu->Diagnostics = malloc(size_diags * sizeof(real)); */

  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    DtFields(simu, simu->w, simu->dtw);
    RK_out(wnp1, simu->w, simu->dtw, 0.5 * dt, simu->wsize);

    simu->tnow += 0.5 * dt;

    DtFields(simu, wnp1, simu->dtw);
    RK_in(simu->w, simu->dtw, dt, simu->wsize);

    simu->tnow += 0.5 * dt;

    /* if(simu->update_after_rk != NULL) */
    /*   simu->update_after_rk(f, simu->wn); */

    iter++;
    simu->iter_time_rk = iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);
  free(wnp1);
}



// Time integration by a fourth-order Runge-Kutta algorithm
void RK4(Simulation *simu, real tmax)
{

  simu->dt = Get_Dt_RK(simu);

  real dt = simu->dt;

  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  // Allocate memory for RK time-stepping
  real *l1, *l2, *l3;
  l1 = calloc(simu->wsize, sizeof(real));
  l2 = calloc(simu->wsize, sizeof(real));
  l3 = calloc(simu->wsize, sizeof(real));
  
  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;
  
    if(simu->nb_diags != 0)
    simu->Diagnostics = malloc(size_diags * sizeof(real));
  
  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    // l_1 = w_n + 0.5dt * S(w_n, t_0)
    DtFields(simu, simu->w, simu->dtw);
    RK_out(l1, simu->w, simu->dtw, 0.5 * dt, simu->wsize);

    simu->tnow += 0.5 * dt;

    // l_2 = w_n + 0.5dt * S(l_1, t_0 + 0.5 * dt)
    DtFields(simu, l1, simu->dtw);
    RK_out(l2, simu->w, simu->dtw, 0.5 * dt, simu->wsize);

    // l_3 = w_n + dt * S(l_2, t_0 + 0.5 * dt)
    DtFields(simu, l2, simu->dtw);
    RK_out(l3, simu->w, simu->dtw, dt, simu->wsize);

    simu->tnow += 0.5 * dt;

    // Compute S(l_3, t_0 + dt)
    DtFields(simu, l3, simu->dtw);
    RK4_final_inplace(simu->w, l1, l2, l3, simu->dtw, dt, simu->wsize);

    
    /* if(f->update_after_rk != NULL) */
    /*   f->update_after_rk(f, f->wn); */
    
    iter++;
     simu->iter_time_rk=iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

  free(l3);
  free(l2);
  free(l1);
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


real Get_Dt_RK(Simulation *simu)
{
  return simu->cfl * simu->hmin / simu->vmax; 
  
}





