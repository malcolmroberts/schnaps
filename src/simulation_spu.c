#include "simulation_spu.h"
#include <stdlib.h>

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

void AddBuffer_C(void *buffers[], void *cl_arg);

void AddBuffer_SPU(real alpha, starpu_data_handle_t win, starpu_data_handle_t wout){

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet AddBuffer...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = AddBuffer_C;
    codelet.nbuffers = 2;
    codelet.modes[0] = STARPU_R;
    codelet.modes[1] = STARPU_RW;
    codelet.name="AddBuffer";
  }

  void* arg_buffer;
  size_t arg_buffer_size;

  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, &alpha, sizeof(real),
			   0);   

  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  task->handles[0] = win;
  task->handles[1] = wout;
    

  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
  
}


void AddBuffer_C(void *buffers[], void *cl_args) {

  real alpha;

  starpu_codelet_unpack_args(cl_args,&alpha);
  free(cl_args);
  
  struct starpu_vector_interface *win_v =
    (struct starpu_vector_interface *) buffers[0];
  real* win = (real *)STARPU_VECTOR_GET_PTR(win_v);

  struct starpu_vector_interface *wout_v =
    (struct starpu_vector_interface *) buffers[1];
  real* wout = (real *)STARPU_VECTOR_GET_PTR(wout_v);

  int n = STARPU_VECTOR_GET_NX(buffers[0]);
  int np = STARPU_VECTOR_GET_NX(buffers[1]);

  assert(n == np);

  for(int i = 0; i < n; i++) wout[i] += alpha * win[i];
  
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
	  f->dtwn[imemL] -= flux[iv]  ;
	  //printf("imem=%d res=%f\n",imemL,f->dtwn[imemL]);
	}
      } else { // The point is on the boundary.


	assert(sign == 1);
	real* xpg = inter->xpg + 3 * ipgf;
	real* vnds = inter->vnds + 3 * ipgf;

	//printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
	f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
	//printf("wL=%f, ipgf=%d\n",wL[0], ipgf);
	//printf("flux=%f, ipgf=%d\n",flux[0], ipgf);
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
    ZeroBuffer_SPU(simu->fd[ie].res_handle);
  }

  starpu_task_wait_for_all();
  
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  }

    // the field pointers must be updated
  if (w_handle != NULL){
    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      simu->fd[ie].wn_handle = w_handle[ie];
      simu->fd[ie].dtwn_handle = dtw_handle[ie];
      simu->fd[ie].res_handle = simu->res_handle[ie];
    }
  }
  else {
    assert(dtw_handle == NULL);
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
	DGMacroCellBoundaryFlux_SPU(inter);
  	//InterfaceBoundaryFlux_SPU(inter);
      }
  }

  
  
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    DGSubCellInterface_SPU(simu->fd + ie);
    DGVolume_SPU(simu->fd + ie);
    DGSource_SPU(simu->fd + ie);
    DGMass_SPU(simu->fd + ie);

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
      rhs[imemL] -= flux[iv] ; ///* f->dt;
      //printf("imem=%d res=%f dt=%f\n",imemL,rhs[imemL],f->dt);
      /* real* xpg = inter->xpg + 3 * ipgf; */
      /* real* vnds = inter->vnds + 3 * ipgf; */
      /* printf("interface flux=%f xpg=%f %f vnds=%f %f\n", */
      /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
    }

  }

  
  


}

void DGMacroCellBoundaryFlux_C(void* buffer[], void* cl_args);

void DGMacroCellBoundaryFlux_SPU(Interface* inter)
{


  field *f;
  int locfa;
  starpu_data_handle_t index;

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGMacroCellBoundaryFlux...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGMacroCellBoundaryFlux_C;
    codelet.nbuffers = 5;
    codelet.modes[0] = STARPU_R;  // vol_index
    codelet.modes[1] = STARPU_R; // vnds
    codelet.modes[2] = STARPU_R; // xpg
    codelet.modes[3] = STARPU_R; // wL
    codelet.modes[4] = STARPU_RW; // rhs / res
    codelet.name="DGMacroCellBoundaryFlux";
  }

  
  f = inter->fL;
  locfa = inter->locfaL;
  index = inter->vol_indexL_handle;

  
  void* arg_buffer;
  size_t arg_buffer_size;

  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, f, sizeof(field),
			   STARPU_VALUE, &locfa, sizeof(int),
			   0);   
    
  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  task->handles[0] = index;
  task->handles[1] = inter->vnds_handle;
  task->handles[2] = inter->xpg_handle;
  task->handles[3] = inter->wL_handle;
  if (f->solver != NULL){
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[4] = sky_spu->rhs_handle;
  }
  else {
    task->handles[4] = f->res_handle;
  }
  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

  
}

void DGMacroCellBoundaryFlux_C(void* buffer[], void* cl_args){


  field f0;
  field *f = &f0;
  int locfa;

  starpu_codelet_unpack_args(cl_args,f,&locfa);
  free(cl_args);

  int m = f->model.m;

  int buf_num=0;
  
  struct starpu_vector_interface *index_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  int* index = (int *)STARPU_VECTOR_GET_PTR(index_v);  

  struct starpu_vector_interface *vnds_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* vnds_buf = (real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);  

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* xpg_buf = (real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);  

  struct starpu_vector_interface *wL_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* wL_buf = (real *)STARPU_VECTOR_GET_PTR(wL_buf_v);  

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* rhs = (real *)STARPU_VECTOR_GET_PTR(rhs_v);  

  int npgf = NPGF(f->deg, f->raf, locfa);
  for(int ipgf = 0; ipgf < npgf; ipgf++) {

 
    real flux[m];
    real wL[m];
    for(int iv = 0; iv < m; iv++) {
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      wL[iv] = wL_buf[imemf];
    }

    real* xpg = xpg_buf + 3 * ipgf;
    real* vnds = vnds_buf + 3 * ipgf;

    //printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
    f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
    //printf("wL=%f, ipgf=%d\n",wL[0], ipgf);
    //printf("flux=%f, ipgf=%d\n",flux[0], ipgf);
    int ipgL = index[ipgf];
    /* printf("xpg=%f %f %f vnds=%f %f %f ipgL=%d \n", */
    /*        xpg[0], xpg[1], xpg[2], */
    /*        vnds[0], vnds[1],vnds[2], ipgL); */
    for(int iv = 0; iv < m; iv++) {
      int ipgL = index[ipgf];
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      rhs[imemL] -= flux[iv]; 
    }
    

  }


}

void DGSubCellInterface_C(void *buffers[], void *cl_arg);

void DGSubCellInterface_SPU(field *f){

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGSubCellInterface...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGSubCellInterface_C;
    codelet.nbuffers = 2;
    codelet.modes[0] = STARPU_R;
    codelet.modes[1] = STARPU_RW;
    codelet.name="DGSubCellInterface";
  }

  void* arg_buffer;
  size_t arg_buffer_size;
   
  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.NumFlux, sizeof(fluxptr),
			   0);


  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;
    

  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
  
}


// Compute inter-subcell fluxes
void DGSubCellInterface_C(void *buffers[], void *cl_arg) 
{

  int m;
  int deg[3];
  int raf[3];
  real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;
  
  starpu_codelet_unpack_args(cl_arg,
			     &m, deg, raf, physnode, &Varindex, &NumFlux);

  free(cl_arg);

  int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};

  struct starpu_vector_interface *w_v =
    (struct starpu_vector_interface *) buffers[0]; 
  real* w = (real *)STARPU_VECTOR_GET_PTR(w_v);  

  struct starpu_vector_interface *res_v =
    (struct starpu_vector_interface *) buffers[1]; 
  real* res = (real *)STARPU_VECTOR_GET_PTR(res_v);  


  
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
		  ref_pg_vol(deg, raf, ipgL, xref, &wpg3, NULL);
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
		  // TO DO change the Varindex signature
		  int imemL = Varindex(deg, raf, m, ipgL, iv); 
		  int imemR = Varindex(deg, raf, m, ipgR, iv);
		  // end TO DO
		  wL[iv] = w[imemL];
		  wR[iv] = w[imemR];
		}
		NumFlux(wL, wR, vnds, flux);

		// subcell ref surface glop weight
		real wpg
		  = wglop(deg[dim1], iL[dim1])
		  * wglop(deg[dim2], iL[dim2]);

		/* printf("vnds %f %f %f flux %f wpg %f\n", */
		/* 	 vnds[0], vnds[1], vnds[2], */
		/* 	 flux[0], wpg); */

		// finally distribute the flux on the two sides
		for(int iv = 0; iv < m; iv++) {
		  // TO DO change the Varindex signature
		  int imemL = Varindex(deg, raf, m, ipgL, iv);
		  int imemR = Varindex(deg, raf, m, ipgR, iv);
		  // end TO DO
		  res[imemL] -= flux[iv] * wpg;
		  res[imemR] += flux[iv] * wpg;
		}

	      }  // face yhat loop
	    } // face xhat loop
	  } // endif internal face
	} // dim loop
      } // subcell icl2 loop
    } // subcell icl1 loop
  } // subcell icl0 loop

}

// Compute the Discontinuous Galerkin volume terms, fast version
void DGVolume_C(void *buffers[], void *cl_arg);

void DGVolume_SPU(field* f)
{
  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGVolume...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGVolume_C;
    codelet.nbuffers = 2;
    codelet.modes[0] = STARPU_R;
    codelet.modes[1] = STARPU_RW;
    codelet.name="DGVolume";
  }

  void* arg_buffer;
  size_t arg_buffer_size;
   
  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.NumFlux, sizeof(fluxptr),
			   0);


  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;
    

  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");





}






void DGVolume_C(void *buffers[], void *cl_arg) 
{


  int m;
  int deg[3];
  int raf[3];
  real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;
  
  starpu_codelet_unpack_args(cl_arg,
			     &m, deg, raf, physnode, &Varindex, &NumFlux);

  free(cl_arg);

  int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};

  const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  struct starpu_vector_interface *w_v =
    (struct starpu_vector_interface *) buffers[0]; 
  real* w = (real *)STARPU_VECTOR_GET_PTR(w_v);  

  struct starpu_vector_interface *res_v =
    (struct starpu_vector_interface *) buffers[1]; 
  real* res = (real *)STARPU_VECTOR_GET_PTR(res_v);  

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

	  ref_pg_vol(deg, raf, offsetL + p, xref, &tomega, NULL);
	  xref0[p] = xref[0];
	  xref1[p] = xref[1];
	  xref2[p] = xref[2];
	  omega[p] = tomega;

	  for(int im = 0; im < m; ++im) {
	    imems[pos++] = Varindex(deg,raf,m, offsetL + p, im);
	  }
	}

	// loop in the "cross" in the three directions
	for(int dim0 = 0; dim0 < 3; dim0++) {
	  //for(int dim0 = 0; dim0 < 2; dim0++) {  // TODO : return to 3d !
	  // point p at which we compute the flux

	  for(int p0 = 0; p0 < npg[0]; p0++) {
	    for(int p1 = 0; p1 < npg[1]; p1++) {
	      for(int p2 = 0; p2 < npg[2]; p2++) {
		real wL[m], flux[m];
		int p[3] = {p0, p1, p2};
		int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
		for(int iv = 0; iv < m; iv++) {
		  ///int imemL = varindex(f_interp_param, ie, ipgL, iv);
		  wL[iv] = w[imems[m * (ipgL - offsetL) + iv]];
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
		  /* ref_pg_vol(interp_param+1,ipgL,xrefL, &wpgL, NULL); */

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

		  NumFlux(wL, wL, dphiL, flux);

		  int ipgR = offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		  for(int iv = 0; iv < m; iv++) {
		    int imemR = Varindex(deg,raf,m, ipgR, iv);
		    int temp = m * (ipgR - offsetL) + iv;  
		    assert(imemR == imems[temp]);
		    res[imems[temp]] += flux[iv] * wpgL;
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

// Apply the source term
void DGSource_C(void *buffers[], void *cl_arg);

void DGSource_SPU(field* f)
{

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGSource...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGSource_C;
    codelet.nbuffers = 2;
    codelet.modes[0] = STARPU_R;
    codelet.modes[1] = STARPU_RW;
    codelet.name="DGSource";
  }

  void* arg_buffer;
  size_t arg_buffer_size;
   
  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.Source, sizeof(sourceptr),
			   STARPU_VALUE, &f->tnow, sizeof(real),
			   0);


  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;
    

  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");




}







void DGSource_C(void *buffers[], void *cl_arg) 
{

  int m;
  int deg[3];
  int raf[3];
  real physnode[20][3];
  varindexptr Varindex;
  sourceptr Source;
  real tnow;
  
  starpu_codelet_unpack_args(cl_arg,
			     &m, deg, raf, physnode, &Varindex, &Source, &tnow);

  free(cl_arg);

  if (Source == NULL) {
    return;
  }


  struct starpu_vector_interface *w_v =
    (struct starpu_vector_interface *) buffers[0]; 
  real* w = (real *)STARPU_VECTOR_GET_PTR(w_v);  

  struct starpu_vector_interface *res_v =
    (struct starpu_vector_interface *) buffers[1]; 
  real* res = (real *)STARPU_VECTOR_GET_PTR(res_v);  

  for(int ipg = 0; ipg < NPG(deg, raf); ipg++) {
    real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(deg, raf, ipg, xpgref, &wpg, NULL);
    Ref2Phy(physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    real det = dot_product(dtau[0], codtau[0]);  //// temp !!!!!!!!!!!!!!!
    real wL[m], source[m];
    for(int iv = 0; iv < m; ++iv){
      int imem = Varindex(deg, raf, m, ipg, iv);
      wL[iv] = w[imem];
    }
      
    Source(xphy, tnow, wL, source);
    // printf("tnow=%f\n",tnow);

    for(int iv = 0; iv < m; ++iv) {
      int imem = Varindex(deg, raf, m, ipg, iv);
      res[imem] += source[iv] * det * wpg;
	
    }
  }
  
}


// Apply division by the mass matrix
void DGMass_C(void *buffers[], void *cl_arg);
void DGMass_OCL(void *buffers[], void *cl_arg);

void DGMass_SPU(field* f)
{

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;
  static struct starpu_opencl_program opencl_program;

  printf("STARPU_MAXIMPLEMENTATIONS=%d\n",STARPU_MAXIMPLEMENTATIONS);
  printf("STARPU_USE_OPENCL=%d\n",STARPU_USE_OPENCL);
  
  if (!is_init){
    printf("init codelet DGMass...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGMass_C;
    codelet.cpu_funcs[1] = NULL;
    printf("Coucou: %s\n",cl_buildoptions);
    int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
    					      &opencl_program, cl_buildoptions);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
    assert(1==3);
    codelet.opencl_funcs[0] =DGMass_OCL;
    codelet.opencl_funcs[1] = NULL;
    codelet.nbuffers = 2;
    codelet.modes[0] = STARPU_R;
    codelet.modes[1] = STARPU_RW;
    codelet.name="DGMass";
  }

  void* arg_buffer;
  size_t arg_buffer_size;
   
  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   0);


  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  task->handles[0] = f->res_handle;
  task->handles[1] = f->dtwn_handle;
    

  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");



}






void DGMass_C(void *buffers[], void *cl_arg) 
{

  int m;
  int deg[3];
  int raf[3];
  real physnode[20][3];
  varindexptr Varindex;
  
  starpu_codelet_unpack_args(cl_arg,
			     &m, deg, raf, physnode, &Varindex);

  free(cl_arg);


  struct starpu_vector_interface *res_v =
    (struct starpu_vector_interface *) buffers[0]; 
  real* res = (real *)STARPU_VECTOR_GET_PTR(res_v);  

  struct starpu_vector_interface *dtw_v =
    (struct starpu_vector_interface *) buffers[1]; 
  real* dtw = (real *)STARPU_VECTOR_GET_PTR(dtw_v);  

  for(int ipg = 0; ipg < NPG(deg, raf); ipg++) {
    real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(deg, raf, ipg, xpgref, &wpg, NULL);
    Ref2Phy(physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    real det = dot_product(dtau[0], codtau[0]);
    for(int iv = 0; iv < m; iv++) {
      int imem = Varindex(deg, raf, m, ipg, iv);
      dtw[imem] = res[imem]/(wpg * det) - dtw[imem] / 2 ;
    }
  }
  
}

void DGMass_OCL(void *buffers[], void *cl_arg) 
{

  int m;
  int deg[3];
  int raf[3];
  real physnode[20][3];
  varindexptr Varindex;
  
  starpu_codelet_unpack_args(cl_arg,
			     &m, deg, raf, physnode, &Varindex);

  free(cl_arg);


  struct starpu_vector_interface *res_v =
    (struct starpu_vector_interface *) buffers[0]; 
  real* res = (real *)STARPU_VECTOR_GET_PTR(res_v);  

  struct starpu_vector_interface *dtw_v =
    (struct starpu_vector_interface *) buffers[1]; 
  real* dtw = (real *)STARPU_VECTOR_GET_PTR(dtw_v);  

  for(int ipg = 0; ipg < NPG(deg, raf); ipg++) {
    real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(deg, raf, ipg, xpgref, &wpg, NULL);
    Ref2Phy(physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    real det = dot_product(dtau[0], codtau[0]);
    for(int iv = 0; iv < m; iv++) {
      int imem = Varindex(deg, raf, m, ipg, iv);
      dtw[imem] = res[imem]/(wpg * det) - dtw[imem] / 2 ;
    }
  }
  
}



void RK2_SPU(Simulation *simu, real tmax){

  simu->dt = Get_Dt_RK(simu);
  real dt = simu->dt;

  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt + 1;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;
  
  // FIXME: remove
  //size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

   /* if(simu->nb_diags != 0) { */
   /*   simu->Diagnostics = malloc(size_diags * sizeof(real)); */
   /* } */

  simu->tnow = 0;

  assert(starpu_use);
  RegisterSimulation_SPU(simu);


  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      ZeroBuffer_SPU(simu->fd[ie].dtwn_handle);
    }
    
    DtFields_SPU(simu, NULL, NULL);
    //UnregisterSimulation_SPU(simu);
    /* for(int i=0 ; i<simu->wsize; i++) */
    /*   printf("i=%d dtw=%f\n",i,simu->dtw[i]); */
    /* assert(1==2); */

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      real alpha = dt/2;
      AddBuffer_SPU(alpha, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle);
    }
    
    //UnregisterSimulation_SPU(simu);
    /* for(int i=0 ; i<simu->wsize; i++) */
    /*   //printf("i=%d w=%f\n",i,simu->w[i]+dt/2*simu->dtw[i]); */
    /*   printf("i=%d w=%f\n",i,simu->w[i]); */
    /* assert(1==2); */
    simu->tnow += 0.5 * dt;

    DtFields_SPU(simu, NULL, NULL);

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      AddBuffer_SPU(simu->dt, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle);
    }

    simu->tnow += 0.5 * dt;

    /* if(simu->update_after_rk != NULL){  */
    /*   simu->update_after_rk(simu, simu->w);  */
    /* } */
    
    iter++;
    simu->iter_time_rk = iter;
  }

  UnregisterSimulation_SPU(simu);

  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);
}
