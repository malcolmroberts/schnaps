#include "simulation_spu.h"


void ZeroBuffer_C(void *buffers[], void *cl_arg);

void ZeroBuffer_SPU(starpu_data_handle_t w){

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet ZeroBuffer...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = ZeroBuffer_C;
    codelet.nbuffers = 1;
    codelet.modes[0] = STARPU_W;
    codelet.name="Zero Buffer";
  }


  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = NULL;
  task->cl_arg_size = 0;
  task->handles[0] = w;
    

  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
  
}


void ZeroBuffer_C(void *buffers[], void *cl_arg) {

  struct starpu_vector_interface *w_v =
    (struct starpu_vector_interface *) buffers[0];

  int n = STARPU_VECTOR_GET_NX(buffers[0]);
  
  real* w = (real *)STARPU_VECTOR_GET_PTR(w_v);

  for(int i = 0; i < n; i++) w[i] = 0;
  
}


void InterfaceExplicitFlux_bis(Interface* inter, int side){


  field *f;
  field *fext;
  int locfa;
  int *index_ext;
  int *index;
  real* wn;
  real* wn_ext;

  const int sign = 1 - 2 * side;
  
  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    index_ext = inter->vol_indexR;
    index = inter->vol_indexL;
    wn = inter->wL;
    wn_ext = inter->wR;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    index_ext = inter->vol_indexL;
    index = inter->vol_indexR;
    wn = inter->wR;
    wn_ext = inter->wL;
  }
  else {
    assert(1==2);
  }

  int npgf = NPGF(f->deg, f->raf, locfa);

  if (f != NULL){
  
    const unsigned int m = f->model.m;

    //printf("locfa=%d \n",locfa);

    for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {

 
      real flux[m];
      real wL[m];
      int ipgL = index[ipgf];
      for(int iv = 0; iv < m; iv++) {
	int imemf = VarindexFace(npgf, m, ipgf, iv);
	wL[iv] = wn[imemf];
	  //wL[iv] = 0;
      }
 
      if (fext != NULL) {  // the right element exists
	real wR[m];
	int ipgR = index_ext[ipgf];
	//int ipgL = index[ipgf];
	for(int iv = 0; iv < m; iv++) {
	  int imemf = VarindexFace(npgf, m, ipgf, iv);
	  wR[iv] = wn_ext[imemf];
	  //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}

	// int_dL F(wL, wR, grad phi_ib)
	real vndsloc[3];
	
	vndsloc[0] = sign * inter->vnds[3 * ipgf + 0];
	vndsloc[1] = sign * inter->vnds[3 * ipgf + 1];
	vndsloc[2] = sign * inter->vnds[3 * ipgf + 2];

	f->model.NumFlux(wL, wR, vndsloc, flux);
 	//printf("flux=%f %f\n",flux[0],flux[0]);

	// Add flux  to the selected side 
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  f->dtwn[imemL] -= flux[iv];
	}
      } else { // The point is on the boundary.


	assert(sign == 1);
	real* xpg = inter->xpg + 3 * ipgf;
	real* vnds = inter->vnds + 3 * ipgf;

	//printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
	f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
	int ipgL = index[ipgf];
	/* printf("xpg=%f %f %f vnds=%f %f %f ipgL=%d \n", */
	/*        xpg[0], xpg[1], xpg[2], */
	/*        vnds[0], vnds[1],vnds[2], ipgL); */
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  f->dtwn[imemL] -= flux[iv] * sign;
	  /* printf("boundary flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
	}
      }

    }

  
  }
}


void DtFields_bis(Simulation *simu,
		  real* w,
		  real* dtw){

  if(simu->pre_dtfields != NULL) {
    simu->pre_dtfields(simu, w);
  }


#ifdef _OPENMP
#pragma omp parallel
#endif

  int fsize =  simu->wsize / simu->macromesh.nbelems;
  
  for(int iw = 0; iw < simu->wsize; iw++)
    dtw[iw] = 0;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  }

    // the field pointers must be updated
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].wn = w + ie * fsize;
    simu->fd[ie].dtwn = dtw + ie * fsize;
  }


  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    // left = 0  right = 1
    ExtractInterface(inter, 0);
    ExtractInterface(inter, 1);
    if (inter->fR != NULL) {
      InterfaceExplicitFlux_bis(inter, 0);
      InterfaceExplicitFlux_bis(inter, 1);
      }
      else{
	InterfaceExplicitFlux_bis(inter, 0);
      //InterfaceBoundaryFlux(inter);
      }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    DGSubCellInterface(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGVolume(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGSource(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGMass(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);

  }

  if(simu->post_dtfields != NULL) {
    simu->post_dtfields(simu, w);
  }

}


void DtFields_SPU(Simulation *simu,
		  starpu_data_handle_t* w_handle,
		  starpu_data_handle_t* dtw_handle)
{

  /* if(simu->pre_dtfields != NULL) { */
  /*   simu->pre_dtfields(simu, w); */
  /* } */

  int fsize =  simu->wsize / simu->macromesh.nbelems;
  
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    ZeroBuffer_SPU(dtw_handle[ie]);
  }

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  }

    // the field pointers must be updated
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].wn_handle = w_handle[ie];
    simu->fd[ie].dtwn_handle = dtw_handle[ie];
    simu->fd[ie].res_handle = dtw_handle[ie];
  }


  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    // left = 0  right = 1
    ExtractInterface_SPU(inter, 0);
    ExtractInterface_SPU(inter, 1);
    if (inter->fR != NULL) {
      InterfaceExplicitFlux_SPU(inter, 0);
      InterfaceExplicitFlux_SPU(inter, 1);
      }
      else{
	InterfaceBoundaryFlux_SPU(inter);
      }
  }

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    /* DGSubCellInterface_SPU(simu->fd + ie); */
    /* DGVolume_SPU(simu->fd + ie); */
    /* DGSource_SPU(simu->fd + ie); */
    /* DGMass_SPU(simu->fd + ie); */

  }

}






