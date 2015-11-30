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
	  printf("imem=%d res=%f\n",imemL,f->dtwn[imemL]);
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
	//InterfaceExplicitFlux_bis(inter, 0);
      //InterfaceBoundaryFlux(inter);
      }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    /* DGSubCellInterface(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */
    /* DGVolume(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */
    /* DGSource(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */
    /* DGMass(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */

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
    simu->fd[ie].res_handle = simu->res_handle[ie];
  }


  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    // left = 0  right = 1
    ExtractInterface_SPU(inter, 0);
    ExtractInterface_SPU(inter, 1);
    if (inter->fR != NULL) {
      DGMacroCellInterface_SPU(inter, 0);
      DGMacroCellInterface_SPU(inter, 1);
      }
      else{
  	//InterfaceBoundaryFlux_SPU(inter);
      }
  }

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    /* DGSubCellInterface_SPU(simu->fd + ie); */
    /* DGVolume_SPU(simu->fd + ie); */
    /* DGSource_SPU(simu->fd + ie); */
    /* DGMass_SPU(simu->fd + ie); */

  }

}


void DGMacroCellInterface_C(void* buffer[], void* cl_args);

void DGMacroCellInterface_SPU(Interface* inter, int side)
{


  field *f;
  field *fext;
  int locfa;
  starpu_data_handle_t index_ext;
  starpu_data_handle_t index;
  starpu_data_handle_t wn_in;
  starpu_data_handle_t wn_ext;

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGMacroCellInterface...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGMacroCellInterface_C;
    codelet.nbuffers = 7;
    codelet.modes[0] = STARPU_R;  // index_ext
    codelet.modes[1] = STARPU_R; // index
    codelet.modes[2] = STARPU_R; // vnds
    codelet.modes[3] = STARPU_R; // xpg
    codelet.modes[4] = STARPU_R; // wn_in
    codelet.modes[5] = STARPU_R; // wn_ext
    codelet.modes[6] = STARPU_RW; // rhs
    codelet.name="DGMacroCellInterface";
  }

  const int sign = 1 - 2 * side;
  
  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    index_ext = inter->vol_indexR_handle;
    wn_in = inter->wL_handle;
    wn_ext = inter->wR_handle;
    index = inter->vol_indexL_handle;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    index_ext = inter->vol_indexL_handle;
    wn_in = inter->wR_handle;
    wn_ext = inter->wL_handle;
    index = inter->vol_indexR_handle;
  }
  else {
    assert(1==2);
  }


  if (f != NULL && fext != NULL){
  
    const unsigned int m = f->model.m;

    void* arg_buffer;
    size_t arg_buffer_size;

    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			     STARPU_VALUE, f, sizeof(field),
			     STARPU_VALUE, fext, sizeof(field),
			     STARPU_VALUE, &locfa, sizeof(int),
			     STARPU_VALUE, &sign, sizeof(int),
			     0);   
    //printf("sizeof_field=%d\n", (int) sizeof(field));
    
    task = starpu_task_create();
    task->cl = &codelet;
    task->cl_arg = arg_buffer;
    task->cl_arg_size = arg_buffer_size;
    task->handles[0] = index_ext;
    task->handles[1] = index;
    task->handles[2] = inter->vnds_handle;
    task->handles[3] = inter->xpg_handle;
    task->handles[4] = wn_in;
    task->handles[5] = wn_ext;
    if (f->solver != NULL){
      Skyline_SPU* sky_spu = f->solver->matrix;
      task->handles[6] = sky_spu->rhs_handle;
    }
    else {
      task->handles[6] = f->res_handle;
    }
    int ret = starpu_task_submit(task);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

  
  }
}


void DGMacroCellInterface_C(void* buffer[], void* cl_args){


  field f0;
  field fext0;
  field *f = &f0;
  field *fext = &fext0;
  int locfa;
  int sign;

  starpu_codelet_unpack_args(cl_args,f,fext,&locfa,&sign);
  free(cl_args);


  int buf_num=0;
  
  struct starpu_vector_interface *index_ext_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  int* index_ext = (int *)STARPU_VECTOR_GET_PTR(index_ext_v);  

  struct starpu_vector_interface *index_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  int* index = (int *)STARPU_VECTOR_GET_PTR(index_v);  

  struct starpu_vector_interface *vnds_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* vnds_buf = (real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);  

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* xpg_buf = (real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);  

  struct starpu_vector_interface *wn_in_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* wn_in = (real *)STARPU_VECTOR_GET_PTR(wn_in_v);  

  struct starpu_vector_interface *wn_ext_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* wn_ext = (real *)STARPU_VECTOR_GET_PTR(wn_ext_v);  

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* rhs = (real *)STARPU_VECTOR_GET_PTR(rhs_v);  
  

  int npgf = NPGF(f->deg, f->raf, locfa);
  int m = fext->model.m;
  
  for(int ipgf = 0; ipgf < npgf; ipgf++) {

 
    real flux[m];
    real wL[m];
    real wR[m];
    int ipgR = index_ext[ipgf];
    int ipgL = index[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      //int imemR = fext->varindex(fext->deg, fext->raf,m , ipgR, iv);
      wL[iv] = wn_in[imemf];
      wR[iv] = wn_ext[imemf];
      //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    // int_dL F(wL, wR, grad phi_ib)
    real vndsloc[3];
	
    vndsloc[0] = sign * vnds_buf[3 * ipgf + 0];
    vndsloc[1] = sign * vnds_buf[3 * ipgf + 1];
    vndsloc[2] = sign * vnds_buf[3 * ipgf + 2];

    //printf("sign=%d ipgL=%d ipgR=%d vndsloc=%f %f\n",sign,ipgL,ipgR,vndsloc[0],vndsloc[1]);

    f->model.NumFlux(wL, wR, vndsloc, flux);
    //printf("flux=%f %f\n",flux[0],flux[0]);

    // Add flux  to the selected side 
    for(int iv = 0; iv < m; iv++) {
      int ipgL = index[ipgf];
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      rhs[imemL] -= flux[iv] * sign; ///* f->dt;
      printf("imem=%d res=%f dt=%f\n",imemL,rhs[imemL],f->dt);
      /* real* xpg = inter->xpg + 3 * ipgf; */
      /* real* vnds = inter->vnds + 3 * ipgf; */
      /* printf("interface flux=%f xpg=%f %f vnds=%f %f\n", */
      /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
    }

  }

  
  


}




