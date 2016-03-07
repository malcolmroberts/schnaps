#include "interface.h"
#include <assert.h>
#include <stdio.h>

#pragma start_opencl
int VarindexFace(int npg, int m, int ipgf, int iv){

  return ipgf * m + iv;

}
#pragma end_opencl

void InitInterface_SPU(Interface* inter){
  if (!starpu_is_init && starpu_use){
    assert(starpu_init(NULL) != -ENODEV);
    starpu_is_init = true;
  }

  if (starpu_use) {
    starpu_vector_data_register(&(inter->vol_indexL_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->vol_indexL), // vector location
				inter->npgL,  // size
				sizeof(int));  // type

    starpu_vector_data_register(&(inter->vol_indexR_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->vol_indexR), // vector location
				inter->npgR,  // size
				sizeof(int));  // type

    starpu_vector_data_register(&(inter->wL_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->wL), // vector location
				inter->wsizeL,  // size
				sizeof(schnaps_real));  // type

    starpu_vector_data_register(&(inter->wR_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->wR), // vector location
				inter->wsizeR,  // size
				sizeof(schnaps_real));  // type

    starpu_vector_data_register(&(inter->vnds_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->vnds), // vector location
				inter->npgL * 3,  // size  !!!!!!!!!!! same for left and right ????
				sizeof(schnaps_real));  // type

    starpu_vector_data_register(&(inter->xpg_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->xpg), // vector location
				inter->npgL * 3,  // size  !!!!!!!!!!! same for left and right ????
				sizeof(schnaps_real));  // type

    starpu_vector_data_register(&(inter->wpg_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->wpg), // vector location
				inter->npgL,  // size  !!!!!!!!!!! same for left and right ????
				sizeof(schnaps_real));  // type

    inter->starpu_registered = true;
  }
}


void UnregisterInterface_SPU(Interface* inter) {
  if (starpu_use && inter->starpu_registered) {
    starpu_data_unregister(inter->vol_indexL_handle);
    starpu_data_unregister(inter->vol_indexR_handle);
    starpu_data_unregister(inter->wL_handle);
    starpu_data_unregister(inter->wR_handle);
    starpu_data_unregister(inter->vnds_handle);
    starpu_data_unregister(inter->xpg_handle);
    starpu_data_unregister(inter->wpg_handle);

    inter->starpu_registered = false;
  }
}


void ExtractInterface_C(void* buffer[], void* cl_args);
void ExtractInterface_OCL(void* buffer[], void* cl_args);

void ExtractInterface_SPU(Interface* inter, int side){

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet ExtractInterface...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = ExtractInterface_C;
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
      codelet.opencl_funcs[0] = ExtractInterface_OCL;
      codelet.opencl_funcs[1] = NULL;
    }
    codelet.nbuffers = 3;
    codelet.modes[0] = STARPU_W;  // wface
    codelet.modes[1] = STARPU_R; // wvol
    codelet.modes[2] = STARPU_R; // vol_index
    codelet.name="ExtractInterface";
  }


    field *fd;
    int locfa,npgf;
    schnaps_real *wf;
    starpu_data_handle_t wf_handle;
    int* vol_index;
    starpu_data_handle_t vol_index_handle;

    if (side == 0) {
      if (inter->fL == NULL) return;
      fd = inter->fL;
      locfa = inter->locfaL;
      npgf = inter->npgL;
      wf = inter->wL;
      wf_handle = inter->wL_handle;
      vol_index = inter->vol_indexL;
      vol_index_handle = inter->vol_indexL_handle;
    }

    if (side == 1) {
      if (inter->fR == NULL) return;
      fd = inter->fR;
      locfa = inter->locfaR;
      npgf = inter->npgR;
      wf = inter->wR;
      wf_handle = inter->wR_handle;
      vol_index = inter->vol_indexR;
      vol_index_handle = inter->vol_indexR_handle;
    }

    assert(side ==1 || side == 0);



    void* arg_buffer;
    size_t arg_buffer_size;

    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			     STARPU_VALUE, &fd->model.m, sizeof(int),
			     STARPU_VALUE, fd->deg, 3 * sizeof(int),
			     STARPU_VALUE, fd->raf, 3 * sizeof(int),
			     STARPU_VALUE, &fd->varindex, sizeof(varindexptr),
			     STARPU_VALUE, &npgf, sizeof(int),
			     0);


    task = starpu_task_create();
    task->cl = &codelet;
    task->cl_arg = arg_buffer;
    task->cl_arg_size = arg_buffer_size;
    task->handles[0] = wf_handle;
    /* Skyline_SPU* sky = fd->solver->matrix; */
    /* task->handles[1] = sky->sol_handle; */
    if (fd->solver != NULL){
      Skyline_SPU* sky_spu = fd->solver->matrix;
      task->handles[1] = sky_spu->sol_handle;
    }
    else {
      task->handles[1] = fd->wn_handle;
    }
    task->handles[2] = vol_index_handle;
    int ret = starpu_task_submit(task);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");




}

void ExtractInterface_C(void* buffer[], void* cl_args){

  int npgf;
  int m;
  int deg[3];
  int raf[3];
  varindexptr Varindex;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, &Varindex,
			     &npgf);
  free(cl_args);


  struct starpu_vector_interface *wf_v =
    (struct starpu_vector_interface *) buffer[0];
  schnaps_real* wf = (schnaps_real *)STARPU_VECTOR_GET_PTR(wf_v);

  struct starpu_vector_interface *wn_v =
    (struct starpu_vector_interface *) buffer[1];
  schnaps_real* wn = (schnaps_real *)STARPU_VECTOR_GET_PTR(wn_v);

  struct starpu_vector_interface *vol_index_v =
    (struct starpu_vector_interface *) buffer[2];
  int* vol_index = (int *)STARPU_VECTOR_GET_PTR(vol_index_v);

  for(int ipgf = 0; ipgf < npgf; ipgf++){
    int ipgv = vol_index[ipgf];
    for(int iv=0; iv < m; iv++){
      int imem = Varindex(deg, raf, m,
			  ipgv, iv);
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      wf[imemf] = wn[imem];
      //printf("wf=%f imem=%d\n",wf[imemf],imem);
    }
  }

}


void ExtractInterface_OCL(void* buffers[], void* cl_args){

  int npgf;
  int m;
  int deg[3];
  int raf[3];
  varindexptr Varindex;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, &Varindex,
			     &npgf);
  free(cl_args);

  int buf_num=0;
  cl_mem wf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem wn = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);
  cl_mem index = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[buf_num++]);


  int id = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(id);
  cl_context context;
  starpu_opencl_get_context(devid, &context);
  cl_event event;

  cl_kernel kernel;
  cl_command_queue queue;

  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program, "ExtractInterface",
                                            devid);
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


  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(int), &m);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &deg_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &raf_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &index);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wf);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
  assert(status >= CL_SUCCESS);

  // Number of glops on the macrocell interface
  size_t wsize = npgf;

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


void ExtractInterface(Interface* inter, int side){

  field *fd;
  int locfa,npgf;
  schnaps_real *wf;
  int mod_ie;
  int* vol_index;

  if (side == 0) {
    fd = inter->fL;
    locfa = inter->locfaL;
    mod_ie = inter->ieL;
    npgf = inter->npgL;
    wf = inter->wL;
    vol_index = inter->vol_indexL;
  }

  if (side == 1) {
    fd = inter->fR;
    locfa = inter->locfaR;
    mod_ie = inter->ieR;
    npgf = inter->npgR;
    wf = inter->wR;
    vol_index = inter->vol_indexR;
  }

  assert(side ==1 || side == 0);

  if (fd !=NULL){
    for(int ipgf = 0; ipgf < npgf; ipgf++){
      int ipgv = vol_index[ipgf];
      for(int iv=0; iv < fd->model.m; iv++){
	int imem = fd->varindex(fd->deg, fd->raf, fd->model.m,
				ipgv, iv);
	int imemf = VarindexFace(npgf, fd->model.m, ipgf, iv);
	//printf("extract to=%d \n",mod_ie);
	wf[imemf] = fd->wn[imem];
      }
    }
  }


}

void InterfaceExplicitFlux_C(void* buffer[], void* cl_args);

void InterfaceExplicitFlux_SPU(Interface* inter, int side)
{


  field *f;
  field *fext;
  int locfa;
  starpu_data_handle_t index_ext;
  starpu_data_handle_t index;
  starpu_data_handle_t wn_ext;

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet InterfaceExplicitFlux...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = InterfaceExplicitFlux_C;
    codelet.nbuffers = 7;
    codelet.modes[0] = STARPU_R;  // index_ext
    codelet.modes[1] = STARPU_R; // index
    codelet.modes[2] = STARPU_R; // vnds
    codelet.modes[3] = STARPU_R; // xpg
    codelet.modes[4] = STARPU_R; // wpg
    codelet.modes[5] = STARPU_R; // wn_ext
    codelet.modes[6] = STARPU_RW; // rhs
    codelet.name="InterfaceExplicitFlux";
  }

  const int sign = 1 - 2 * side;

  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    index_ext = inter->vol_indexR_handle;
    wn_ext = inter->wR_handle;
    index = inter->vol_indexL_handle;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    index_ext = inter->vol_indexL_handle;
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
    int nhandle = 0;
    task->handles[nhandle++] = index_ext;
    task->handles[nhandle++] = index;
    task->handles[nhandle++] = inter->vnds_handle;
    task->handles[nhandle++] = inter->xpg_handle;
    task->handles[nhandle++] = inter->wpg_handle;
    task->handles[nhandle++] = wn_ext;
    if (f->solver != NULL){
      Skyline_SPU* sky_spu = f->solver->matrix;
      task->handles[nhandle++] = sky_spu->rhs_handle;
    }
    else {
      task->handles[nhandle++] = f->res_handle;
    }
    int ret = starpu_task_submit(task);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");


  }
}


void InterfaceExplicitFlux_C(void* buffer[], void* cl_args){


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
  schnaps_real* vnds_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  schnaps_real* xpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);

  struct starpu_vector_interface *wpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  schnaps_real* wpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(wpg_buf_v);

  struct starpu_vector_interface *wn_ext_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  schnaps_real* wn_ext = (schnaps_real *)STARPU_VECTOR_GET_PTR(wn_ext_v);

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  schnaps_real* rhs = (schnaps_real *)STARPU_VECTOR_GET_PTR(rhs_v);


  int npgf = NPGF(f->deg, f->raf, locfa);
  int m = fext->model.m;

  for(int ipgf = 0; ipgf < npgf; ipgf++) {


    schnaps_real flux[m];
    schnaps_real wL[m];
    for(int iv = 0; iv < m; iv++) {
      wL[iv] = 0;
    }

    schnaps_real wR[m];
    int ipgR = index_ext[ipgf];
    int ipgL = index[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      //int imemR = fext->varindex(fext->deg, fext->raf,m , ipgR, iv);
      wR[iv] = wn_ext[imemf];
      //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    // int_dL F(wL, wR, grad phi_ib)
    schnaps_real vndsloc[3];

    vndsloc[0] = sign * vnds_buf[3 * ipgf + 0];
    vndsloc[1] = sign * vnds_buf[3 * ipgf + 1];
    vndsloc[2] = sign * vnds_buf[3 * ipgf + 2];

    //printf("sign=%d ipgL=%d ipgR=%d vndsloc=%f %f\n",sign,ipgL,ipgR,vndsloc[0],vndsloc[1]);

    f->model.NumFlux(wL, wR, vndsloc, flux);
    //printf("flux=%f %f\n",flux[0],flux[0]);

    // Add flux  to the selected side
    const schnaps_real wpg = wpg_buf[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int ipgL = index[ipgf];
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      rhs[imemL] -= flux[iv] * f->dt * wpg;
      printf("imem=%d res=%f\n",imemL,rhs[imemL]);
      /* real* xpg = inter->xpg + 3 * ipgf; */
      /* real* vnds = inter->vnds + 3 * ipgf; */
      /* printf("interface flux=%f xpg=%f %f vnds=%f %f\n", */
      /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
    }

  }





}


void InterfaceBoundaryFlux_C(void* buffer[], void* cl_args);

void InterfaceBoundaryFlux_SPU(Interface* inter)
{


  field *f;
  int locfa;
  starpu_data_handle_t index;

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet InterfaceBoundaryFlux...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = InterfaceBoundaryFlux_C;
    codelet.nbuffers = 5;
    codelet.modes[0] = STARPU_R;  // wface
    codelet.modes[1] = STARPU_R; // wvol
    codelet.modes[2] = STARPU_R; // vol_index
    codelet.modes[3] = STARPU_R; // wpg
    codelet.modes[4] = STARPU_RW; // vol_index
    codelet.name="InterfaceBoundaryFlux";
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
  int nhandle = 0;
  task->handles[nhandle++] = index;
  task->handles[nhandle++] = inter->vnds_handle;
  task->handles[nhandle++] = inter->xpg_handle;
  task->handles[nhandle++] = inter->wpg_handle;
  if (f->solver != NULL){
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[nhandle++] = sky_spu->rhs_handle;
  }
  else {
    task->handles[nhandle++] = f->res_handle;
  }
  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");


}

void InterfaceBoundaryFlux_C(void* buffer[], void* cl_args){


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
  schnaps_real* vnds_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  schnaps_real* xpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);

  struct starpu_vector_interface *wpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  schnaps_real* wpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(wpg_buf_v);

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffer[buf_num++];
  schnaps_real* rhs = (schnaps_real *)STARPU_VECTOR_GET_PTR(rhs_v);

  for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {


    schnaps_real flux[m];
    schnaps_real wL[m];
    for(int iv = 0; iv < m; iv++) {
      wL[iv] = 0;
    }

    schnaps_real* xpg = xpg_buf + 3 * ipgf;
    schnaps_real* vnds = vnds_buf + 3 * ipgf;

    //printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
    f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
    int ipgL = index[ipgf];
    /* printf("xpg=%f %f %f vnds=%f %f %f ipgL=%d \n", */
    /*        xpg[0], xpg[1], xpg[2], */
    /*        vnds[0], vnds[1],vnds[2], ipgL); */
    const schnaps_real wpg = wpg_buf[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int ipgL = index[ipgf];
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      rhs[imemL] -= flux[iv] * f->dt * wpg;
    }


  }








}


void InterfaceExplicitFlux(Interface* inter, int side){


  field *f;
  field *fext;
  int locfa;
  int mod_ie;
  int *index_ext;
  int *index;

  const int sign = 1 - 2 * side;

  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    mod_ie = inter->ieL;
    index_ext = inter->vol_indexR;
    index = inter->vol_indexL;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    mod_ie = inter->ieR;
    index_ext = inter->vol_indexL;
    index = inter->vol_indexR;
  }
  else {
    assert(1==2);
  }


  if (f != NULL){
    schnaps_real* res;
    if (f->solver != NULL){
      res = f->solver->rhs;
    } else {
      res = f->res;
    }

    const unsigned int m = f->model.m;

    //printf("locfa=%d \n",locfa);

    for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {


      schnaps_real flux[m];
      schnaps_real wL[m];
      for(int iv = 0; iv < m; iv++) {
	wL[iv] = 0;
      }

      const schnaps_real wpg = inter->wpg[ipgf];

      if (fext != NULL) {  // the right element exists
	schnaps_real wR[m];
	int ipgR = index_ext[ipgf];
	int ipgL = index[ipgf];
	for(int iv = 0; iv < m; iv++) {
	  int imemR = fext->varindex(fext->deg, fext->raf,fext->model.m, ipgR, iv);
	  wR[iv] = fext->wn[imemR];
	  //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}

	// int_dL F(wL, wR, grad phi_ib)
	schnaps_real vndsloc[3];

	vndsloc[0] = sign * inter->vnds[3 * ipgf + 0];
	vndsloc[1] = sign * inter->vnds[3 * ipgf + 1];
	vndsloc[2] = sign * inter->vnds[3 * ipgf + 2];

	//printf("sign=%d ipgL=%d ipgR=%d vndsloc=%f %f\n",sign,ipgL,ipgR,vndsloc[0],vndsloc[1]);

	//assert(f != NULL);

	f->model.NumFlux(wL, wR, vndsloc, flux);
	//printf("flux=%f %f\n",flux[0],flux[0]);

	// Add flux  to the selected side
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  res[imemL] -= flux[iv] * f->dt * wpg;
	  printf("imem=%d res=%f\n",imemL,res[imemL]);

	  /* real* xpg = inter->xpg + 3 * ipgf; */
	  /* real* vnds = inter->vnds + 3 * ipgf; */
	  /* printf("interface flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
	}
      } else { // The point is on the boundary.


	assert(sign == 1);
	schnaps_real* xpg = inter->xpg + 3 * ipgf;
	schnaps_real* vnds = inter->vnds + 3 * ipgf;

	//printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
	f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
	int ipgL = index[ipgf];
	printf("xpg=%f %f %f vnds=%f %f %f ipgL=%d &rhs=%p\n",
	       xpg[0], xpg[1], xpg[2],
	       vnds[0], vnds[1],vnds[2], ipgL,res);
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  res[imemL] -= flux[iv] * sign * f->dt * wpg;
	  /* printf("boundary flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
	}
      }

    }


  }
}
