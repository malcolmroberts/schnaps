#include "simulation_spu.h"
#include <stdlib.h>


void DisplayHandle_SPU(starpu_data_handle_t handle,
                       const char* name) {
  starpu_task_wait_for_all();
  // Force update of memory on main node
  starpu_data_prefetch_on_node(handle, 0, 0);
  // Warning: only real arrays for the moment
  // Extention to every type could be made with starpu_vector_get_elemsize
  real* ptr = (real*) starpu_vector_get_local_ptr(handle);
  int size = starpu_vector_get_nx(handle);
  for (int i = 0; i < size; ++i)
    printf("%s[%d]: %f\n", name, i , ptr[i]);
}

void ZeroBuffer_C(void *buffers[], void *cl_arg);
void ZeroBuffer_OCL(void *buffers[], void *cl_arg);

void ZeroBuffer_SPU(starpu_data_handle_t w){

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet ZeroBuffer...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    if (starpu_c_use) {
      codelet.cpu_funcs[0] = ZeroBuffer_C;
      codelet.cpu_funcs[1] = NULL;
    }
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = ZeroBuffer_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
    codelet.nbuffers = 1;
    codelet.modes[0] = STARPU_W;
    codelet.name="ZeroBuffer";
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

void ZeroBuffer_OCL(void *buffers[], void *cl_arg) {

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  size_t n = STARPU_VECTOR_GET_NX(buffers[0]);


  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "set_buffer_to_zero",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);

  status = clEnqueueNDRangeKernel(queue,
				  kernel,
				  1, // cl_uint work_dim,
				  NULL, // global_work_offset,
				  &n, // global_work_size,
				  NULL, // local_work_size,
				  0,
				  NULL, // *wait_list,
				  &event); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);
}


void AddBuffer_C(void *buffers[], void *cl_arg);
void AddBuffer_OCL(void *buffers[], void *cl_arg);

// compute wout = wout + alpha * win
void AddBuffer_SPU(real alpha, starpu_data_handle_t win, starpu_data_handle_t wout){

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet AddBuffer...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    if (starpu_c_use) {
      codelet.cpu_funcs[0] = AddBuffer_C;
      codelet.cpu_funcs[1] = NULL;
    }
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = AddBuffer_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
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


void AddBuffer_OCL(void *buffers[], void *cl_args) {

  real alpha;

  starpu_codelet_unpack_args(cl_args,&alpha);
  free(cl_args);

  cl_mem win = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  size_t n = STARPU_VECTOR_GET_NX(buffers[0]);

  cl_mem wout = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);
  size_t np = STARPU_VECTOR_GET_NX(buffers[1]);

  assert(n == np);


  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "AddBuffer",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(real), &alpha);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &win);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wout);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);

  status = clEnqueueNDRangeKernel(queue,
				  kernel,
				  1, // cl_uint work_dim,
				  NULL, // global_work_offset,
				  &n, // global_work_size,
				  NULL, // local_work_size,
				  0,
				  NULL, // *wait_list,
				  &event); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);

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

  //starpu_task_wait_for_all();


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
      // FIXME: significant increase of error when coupling L and R (10**3)
      // Apply to L alone is ok, apply to R alone is ok, apply on L + R is KO
      /* DGMacroCellInterface_bis_SPU(inter); */
      // FIXME: significant increase of error when coupling side 0 and 1 (10**2)
      // Side 0 alone is ok, side 1 alone is ok, side 0 + 1 is KO
      DGMacroCellInterface_SPU(inter, 0);
      DGMacroCellInterface_SPU(inter, 1);
    } else{
      DGMacroCellBoundaryFlux_SPU(inter);
    }
  }


  //printf("STEP ---------------------------------------------------\n");
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    //DisplayHandle_SPU(simu->fd[ie].wn_handle, "wn_INI");
    //DisplayHandle_SPU(simu->fd[ie].res_handle, "res_INI");

    DGSubCellInterface_SPU(simu->fd + ie);

    //DisplayHandle_SPU(simu->fd[ie].res_handle, "res_SCI");

    DGVolume_SPU(simu->fd + ie);

    //DisplayHandle_SPU(simu->fd[ie].res_handle, "res_VOL");

    DGSource_SPU(simu->fd + ie);

    //DisplayHandle_SPU(simu->fd[ie].res_handle, "res_SOU");

    DGMass_SPU(simu->fd + ie);

    //DisplayHandle_SPU(simu->fd[ie].wn_handle, "wn_MAS");

  }

}


void DGMacroCellInterface_bis_C(void* buffer[], void* cl_args);
//void DGMacroCellInterface_OCL(void* buffer[], void* cl_args);

void DGMacroCellInterface_bis_SPU(Interface* inter) {

  field* fL = inter->fL;
  field* fR = inter->fR;

  if (fL != NULL && fR != NULL) {

    static bool is_init = false;
    static struct starpu_codelet codelet;
    struct starpu_task *task;

    if (!is_init) {
      printf("init codelet DGMacroCellInterface_bis...\n");
      is_init = true;
      starpu_codelet_init(&codelet);
      codelet.cpu_funcs[0] = DGMacroCellInterface_bis_C;
      codelet.cpu_funcs[1] = NULL;
      // No kernel "bis" yet, do not uncomment !!!!!!!!!!!!!!!!!!!!!!!
      /* if (starpu_ocl_use) { */
      /*   if (!opencl_program_is_init) { */
      /*     opencl_program_is_init = true; */
      /*     printf("load OpenCL program...\n"); */
      /*     int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl", */
      /*                                                   &opencl_program, */
      /*                                                   cl_buildoptions); */
      /*     STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file"); */
      /*   } */
      /*   codelet.opencl_funcs[0] = DGMacroCellInterface_bis_OCL; */
      /*   codelet.opencl_funcs[1] = NULL; */
      /* } */
      codelet.nbuffers = 7;
      codelet.modes[0] = STARPU_R; // index L
      codelet.modes[1] = STARPU_R; // index R
      codelet.modes[2] = STARPU_R; // vnds
      codelet.modes[3] = STARPU_R; // wn L
      codelet.modes[4] = STARPU_R; // wn R
      codelet.modes[5] = STARPU_RW; // rhs L
      codelet.modes[6] = STARPU_RW; // rhs R
      codelet.name="DGMacroCellInterface_bis";
    }


    void* arg_buffer;
    size_t arg_buffer_size;

    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			     STARPU_VALUE, fL, sizeof(field),
			     STARPU_VALUE, fR, sizeof(field),
                             // We assume that opposite faces have the same npgf!
			     STARPU_VALUE, &(inter->locfaL), sizeof(int),
			     0);

    task = starpu_task_create();
    task->cl = &codelet;
    task->cl_arg = arg_buffer;
    task->cl_arg_size = arg_buffer_size;
    task->handles[0] = inter->vol_indexL_handle;
    task->handles[1] = inter->vol_indexR_handle;
    task->handles[2] = inter->vnds_handle;
    task->handles[3] = inter->wL_handle;
    task->handles[4] = inter->wR_handle;
    if (fL->solver != NULL) {
      Skyline_SPU* sky_spu = fL->solver->matrix;
      task->handles[5] = sky_spu->rhs_handle;
    } else {
      task->handles[5] = fL->res_handle;
    }
    if (fR->solver != NULL) {
      Skyline_SPU* sky_spu = fR->solver->matrix;
      task->handles[6] = sky_spu->rhs_handle;
    } else {
      task->handles[6] = fR->res_handle;
    }
    int ret = starpu_task_submit(task);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
  }
}


void DGMacroCellInterface_bis_C(void* buffer[], void* cl_args) {

  field fL0;
  field fR0;
  field *fL = &fL0;
  field *fR = &fR0;
  int locfa;

  starpu_codelet_unpack_args(cl_args, fL, fR, &locfa);
  free(cl_args);


  int buf_num = 0;

  struct starpu_vector_interface *vol_indexL_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  int* vol_indexL = (int *)STARPU_VECTOR_GET_PTR(vol_indexL_v);

  struct starpu_vector_interface *vol_indexR_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  int* vol_indexR = (int *)STARPU_VECTOR_GET_PTR(vol_indexR_v);

  struct starpu_vector_interface *vnds_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  real* vnds_buf = (real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);

  struct starpu_vector_interface *wnL_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  real* wnL = (real *)STARPU_VECTOR_GET_PTR(wnL_v);

  struct starpu_vector_interface *wnR_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  real* wnR = (real *)STARPU_VECTOR_GET_PTR(wnR_v);

  struct starpu_vector_interface *resL_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  real* resL = (real *)STARPU_VECTOR_GET_PTR(resL_v);

  struct starpu_vector_interface *resR_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  real* resR = (real *)STARPU_VECTOR_GET_PTR(resR_v);


  const int npgf = NPGF(fL->deg, fL->raf, locfa);
  const int m = fL->model.m;

  for (int ipgf = 0; ipgf < npgf; ++ipgf) {

    real wL[m];
    real wR[m];
    for (int iv = 0; iv < m; ++iv) {
      const int imem = VarindexFace(npgf, m, ipgf, iv);
      wL[iv] = wnL[imem];
      wR[iv] = wnR[imem];
    }

    real vnds[3];
    vnds[0] = vnds_buf[3 * ipgf + 0];
    vnds[1] = vnds_buf[3 * ipgf + 1];
    vnds[2] = vnds_buf[3 * ipgf + 2];

    real flux[m];
    fL->model.NumFlux(wL, wR, vnds, flux);

    const int ipgL = vol_indexL[ipgf];
    const int ipgR = vol_indexR[ipgf];
    for (int iv = 0; iv < m; ++iv) {
      const int imemL = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv);
      const int imemR = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, iv);
      resL[imemL] -= flux[iv];
      resR[imemR] += flux[iv];
    }
    /* printf("ipgf=%d ipgL=%d ipgR=%d wL=%f wR=%f vnds=%f %f %f flux=%f\n", */
    /*        ipgf, ipgL, ipgR, wL[0], wR[0], vnds[0], vnds[1], vnds[2], flux[0]); */

  }

}


void DGMacroCellInterface_C(void* buffer[], void* cl_args);
void DGMacroCellInterface_OCL(void* buffer[], void* cl_args);

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
    codelet.cpu_funcs[1] = NULL;
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = DGMacroCellInterface_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
    codelet.nbuffers = 5;
    codelet.modes[0] = STARPU_R; // index
    codelet.modes[1] = STARPU_R; // vnds
    codelet.modes[2] = STARPU_R; // wn_in
    codelet.modes[3] = STARPU_R; // wn_ext
    codelet.modes[4] = STARPU_RW; // rhs
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
    task->handles[0] = index;
    task->handles[1] = inter->vnds_handle;
    task->handles[2] = wn_in;
    task->handles[3] = wn_ext;
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

  struct starpu_vector_interface *index_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  int* index = (int *)STARPU_VECTOR_GET_PTR(index_v);

  struct starpu_vector_interface *vnds_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  real* vnds_buf = (real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);

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

    f->model.NumFlux(wL, wR, vndsloc, flux);
    //printf("flux=%f %f\n",flux[0],flux[0]);

    // Add flux  to the selected side
    int ipgL = index[ipgf];
    for(int iv = 0; iv < m; iv++) {
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      rhs[imemL] -= flux[iv] ; ///* f->dt;
    }
    /* printf("ipgf=%d ipgL=%d wL=%f wR=%f vnds=%f %f %f flux=%f\n", */
    /*        ipgf, ipgL, wL[0], wR[0], vndsloc[0], vndsloc[1], vndsloc[2], flux[0]); */

  }

}

void DGMacroCellInterface_OCL(void* buffers[], void* cl_args){

  field f0;
  field fext0;
  field *f = &f0;
  field *fext = &fext0;
  int locfa;
  int sign;

  starpu_codelet_unpack_args(cl_args, f, fext, &locfa, &sign);
  free(cl_args);

  int buf_num=0;
  cl_mem index = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem vnds_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem wn_in = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem wn_ext = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);


  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "DGMacroCellInterfaceRes",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem deg_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 f->deg,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem raf_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 f->raf,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);


  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(int), &(f->model.m));
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &deg_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &raf_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &sign);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &index);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn_in);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn_ext);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &vnds_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);

  // Number of glops on the macrocell interface
  size_t wsize = NPGF(f->deg, f->raf, locfa);

  status = clEnqueueNDRangeKernel(queue,
                                  kernel,
                                  1, // cl_uint work_dim,
                                  NULL, // global_work_offset,
                                  &wsize, // global_work_size,
                                  NULL, // local_work_size,
                                  0,
                                  NULL, // *wait_list,
                                  &event); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);
}



void DGMacroCellBoundaryFlux_C(void* buffer[], void* cl_args);
void DGMacroCellBoundaryFlux_OCL(void* buffer[], void* cl_args);

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
    codelet.cpu_funcs[1] = NULL;
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = DGMacroCellBoundaryFlux_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
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

void DGMacroCellBoundaryFlux_OCL(void* buffers[], void* cl_args){

  field f0;
  field *f = &f0;
  int locfa;

  starpu_codelet_unpack_args(cl_args,f,&locfa);
  free(cl_args);

  int buf_num=0;
  cl_mem index = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem vnds_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem xpg_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem wn_in = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);

  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "DGBoundaryRes",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem deg_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 f->deg,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem raf_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 f->raf,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);


  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(int), &(f->model.m));
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &deg_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &raf_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(real), &(f->tnow));
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &index);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn_in);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &xpg_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &vnds_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);

  // Number of glops on the macrocell interface
  size_t wsize = NPGF(f->deg, f->raf, locfa);

  status = clEnqueueNDRangeKernel(queue,
                                  kernel,
                                  1, // cl_uint work_dim,
                                  NULL, // global_work_offset,
                                  &wsize, // global_work_size,
                                  NULL, // local_work_size,
                                  0,
                                  NULL, // *wait_list,
                                  &event); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);
}



void DGSubCellInterface_C(void *buffers[], void *cl_arg);
void DGSubCellInterface_OCL(void *buffers[], void *cl_arg);

void DGSubCellInterface_SPU(field *f){

  if (f->raf[0] * f->raf[1] * f->raf[2] == 1)
    return;


  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGSubCellInterface...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGSubCellInterface_C;
    codelet.cpu_funcs[1] = NULL;
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = DGSubCellInterface_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
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

		/* printf("wL %f wR %f vnds %f %f %f flux %f wpg %f\n", */
                /*        wL[0], wR[0], */
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


void DGSubCellInterface_OCL(void *buffers[], void *cl_arg)
{

  int m;
  int deg[3];
  int raf[3];
  real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_arg,
                             &m, deg, raf, physnode,
                             &Varindex, &NumFlux);
  free(cl_arg);

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);


  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "DGSubCellInterfaceRes",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
                                      CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                      sizeof(real) * 60,
                                      physnode,
                                      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem deg_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 deg,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem raf_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 raf,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);


  // Loop over the dimensions
  for (int dim0 = 0; dim0 < 3; ++dim0) {
    // More than one subcell according to that dimension
    if (raf[dim0] > 1) {
      // Other dimensions
      int dim1 = (dim0 + 1) % 3;
      int dim2 = (dim0 + 2) % 3;

      // Number of glops on a subface normal to that dimension (both sides)
      size_t lsize = 2 * (deg[dim1] + 1) * (deg[dim2] + 1);
      // Number subfaces in the macrocell
      size_t wsize = lsize * (raf[dim0] - 1) * raf[dim1] * raf[dim2];

      int narg = 0;
      status = clSetKernelArg(kernel, narg++, sizeof(int), &m);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &deg_cl);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &raf_cl);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
      status |= clSetKernelArg(kernel, narg++, sizeof(int), &dim0);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
      status |= clSetKernelArg(kernel, narg++, sizeof(real) * lsize * m, NULL);
      if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
      assert(status >= CL_SUCCESS);


      status = clEnqueueNDRangeKernel(queue,
                                      kernel,
                                      1, // cl_uint work_dim,
                                      NULL, // global_work_offset,
                                      &wsize, // global_work_size,
                                      &lsize, // local_work_size,
                                      0,
                                      NULL, // *wait_list,
                                      &event); // *event
      if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
      assert(status >= CL_SUCCESS);

    }
  }
}



// Compute the Discontinuous Galerkin volume terms, fast version
void DGVolume_C(void *buffers[], void *cl_arg);
void DGVolume_OCL(void *buffers[], void *cl_arg);

void DGVolume_SPU(field* f)
{
  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGVolume...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    if (starpu_c_use) {
      codelet.cpu_funcs[0] = DGVolume_C;
      codelet.cpu_funcs[1] = NULL;
    }
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = DGVolume_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
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



void DGVolume_OCL(void *buffers[], void *cl_arg) {

  int m;
  int deg[3];
  int raf[3];
  real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_arg,
			     &m, deg, raf, physnode, &Varindex, &NumFlux);

  free(cl_arg);

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);


  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "DGVolumeRes",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
                                      CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                      sizeof(real) * 60,
                                      physnode,
                                      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
                                   CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                   sizeof(int) * 7,
                                   param,
                                   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int ie = 0;

  size_t wsize = NPG(deg, raf);
  size_t lsize = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1);

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  status |= clSetKernelArg(kernel, narg++, sizeof(real) * 2 * lsize * m, NULL);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);


  status = clEnqueueNDRangeKernel(queue,
				  kernel,
				  1, // cl_uint work_dim,
				  NULL, // global_work_offset,
				  &wsize, // global_work_size,
				  &lsize, // local_work_size,
				  0,
				  NULL, // *wait_list,
				  &event); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);


}



// Apply the source term
void DGSource_C(void *buffers[], void *cl_arg);
void DGSource_OCL(void *buffers[], void *cl_arg);

void DGSource_SPU(field* f)
{

  if (f->model.Source == NULL)
    return;


  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGSource...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGSource_C;
    codelet.cpu_funcs[1] = NULL;
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = DGSource_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
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



void DGSource_OCL(void *buffers[], void *cl_arg)
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

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);


  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "DGSourceRes",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
                                      CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                      sizeof(real) * 60,
                                      physnode,
                                      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem deg_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 deg,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem raf_cl = clCreateBuffer(context,
                                 CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                 sizeof(int) * 3,
                                 raf,
                                 &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  size_t wsize = NPG(deg, raf);
  size_t lsize = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1);

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(int), &m);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &deg_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &raf_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(real), &tnow);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  status |= clSetKernelArg(kernel, narg++, sizeof(real) * 2 * lsize * m, NULL);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);


  status = clEnqueueNDRangeKernel(queue,
				  kernel,
				  1, // cl_uint work_dim,
				  NULL, // global_work_offset,
				  &wsize, // global_work_size,
				  &lsize, // local_work_size,
				  0,
				  NULL, // *wait_list,
				  &event); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);


}




// Apply division by the mass matrix
void DGMass_C(void *buffers[], void *cl_arg);
void DGMass_OCL(void *buffers[], void *cl_arg);

void DGMass_SPU(field* f)
{

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet DGMass...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = DGMass_C;
    codelet.cpu_funcs[1] = NULL;
    if (starpu_ocl_use) {
      if (!opencl_program_is_init) {
        opencl_program_is_init = true;
        printf("load OpenCL program...\n");
        int ret = starpu_opencl_load_opencl_from_file("./schnaps.cl",
                                                      &opencl_program,
                                                      cl_buildoptions);
        STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_load_opencl_from_file");
      }
      codelet.opencl_funcs[0] = DGMass_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
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
      dtw[imem] = res[imem] / (wpg * det) - dtw[imem] / 2 ;
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


  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem dtw = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);

  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "DGMassRes",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
			       CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			       sizeof(real) * 60,
			       physnode,
			       &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);


  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(int), &m);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), deg + 0);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), deg + 1);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), deg + 2);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), raf + 0);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), raf + 1);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), raf + 2);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &dtw);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);

  size_t wsize = NPG(deg, raf);

  status = clEnqueueNDRangeKernel(queue,
				  kernel,
				  1, // cl_uint work_dim,
				  NULL, // global_work_offset,
				  &wsize, // global_work_size,
				  NULL, // local_work_size,
				  0,
				  NULL, // *wait_list,
				  &event); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);


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
