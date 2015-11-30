#include "interface.h"
#include <assert.h>
#include <stdio.h>

int VarindexFace(int npg, int m, int ipgf, int iv){

  return ipgf * m + iv;

}

void InitInterface_SPU(Interface* inter){

  if (starpu_use){


  
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
				sizeof(real));  // type

    starpu_vector_data_register(&(inter->wR_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->wR), // vector location
				inter->wsizeR,  // size
				sizeof(real));  // type

    starpu_vector_data_register(&(inter->vnds_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->vnds), // vector location
				inter->npgL * 3,  // size  !!!!!!!!!!! same for left and right ????
				sizeof(real));  // type

    starpu_vector_data_register(&(inter->xpg_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(inter->xpg), // vector location
				inter->npgL * 3,  // size  !!!!!!!!!!! same for left and right ????
				sizeof(real));  // type
  }
}

void ExtractInterface_C(void* buffer[], void* cl_args);
  
void ExtractInterface_SPU(Interface* inter, int side){

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet ExtractInterface...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = ExtractInterface_C;
    codelet.nbuffers = 3;
    codelet.modes[0] = STARPU_W;  // wface
    codelet.modes[1] = STARPU_R; // wvol
    codelet.modes[2] = STARPU_R; // vol_index
    codelet.name="ExtractInterface";
  }


    field *fd;
    int locfa,npgf;
    real *wf;
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
  real* wf = (real *)STARPU_VECTOR_GET_PTR(wf_v);  

  struct starpu_vector_interface *wn_v =
    (struct starpu_vector_interface *) buffer[1]; 
  real* wn = (real *)STARPU_VECTOR_GET_PTR(wn_v);  
  
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


void ExtractInterface(Interface* inter, int side){

  field *fd;
  int locfa,npgf;
  real *wf;
  int* vol_index;

  if (side == 0) {
    fd = inter->fL;
    locfa = inter->locfaL;
    npgf = inter->npgL;
    wf = inter->wL;
    vol_index = inter->vol_indexL;
  }

  if (side == 1) {
    fd = inter->fR;
    locfa = inter->locfaR;
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
    codelet.nbuffers = 6;
    codelet.modes[0] = STARPU_R;  // index_ext
    codelet.modes[1] = STARPU_R; // index
    codelet.modes[2] = STARPU_R; // vnds
    codelet.modes[3] = STARPU_R; // xpg
    codelet.modes[4] = STARPU_R; // wn_ext
    codelet.modes[5] = STARPU_RW; // rhs
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
    task->handles[0] = index_ext;
    task->handles[1] = index;
    task->handles[2] = inter->vnds_handle;
    task->handles[3] = inter->xpg_handle;
    task->handles[4] = wn_ext;
    if (f->solver != NULL){
      Skyline_SPU* sky_spu = f->solver->matrix;
      task->handles[5] = sky_spu->rhs_handle;
    }
    else {
      task->handles[5] = f->res_handle;
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
  real* vnds_buf = (real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);  

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* xpg_buf = (real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);  

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
    for(int iv = 0; iv < m; iv++) {
      wL[iv] = 0;
    }
 
    real wR[m];
    int ipgR = index_ext[ipgf];
    int ipgL = index[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      //int imemR = fext->varindex(fext->deg, fext->raf,m , ipgR, iv);
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
      rhs[imemL] -= flux[iv] * f->dt;
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
    codelet.nbuffers = 4;
    codelet.modes[0] = STARPU_R;  // wface
    codelet.modes[1] = STARPU_R; // wvol
    codelet.modes[2] = STARPU_R; // vol_index
    codelet.modes[3] = STARPU_RW; // vol_index
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
  task->handles[0] = index;
  task->handles[1] = inter->vnds_handle;
  task->handles[2] = inter->xpg_handle;
  if (f->solver != NULL){
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[3] = sky_spu->rhs_handle;
  }
  else {
    task->handles[3] = f->res_handle;
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
  real* vnds_buf = (real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);  

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* xpg_buf = (real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);  

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffer[buf_num++]; 
  real* rhs = (real *)STARPU_VECTOR_GET_PTR(rhs_v);  
  
  for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {

 
    real flux[m];
    real wL[m];
    for(int iv = 0; iv < m; iv++) {
      wL[iv] = 0;
    }
 
    real* xpg = xpg_buf + 3 * ipgf;
    real* vnds = vnds_buf + 3 * ipgf;

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
      rhs[imemL] -= flux[iv] * f->dt;
    }
    

  }

  
  



  

}


void InterfaceExplicitFlux(Interface* inter, int side){


  field *f;
  field *fext;
  int locfa;
  int *index_ext;
  int *index;

  const int sign = 1 - 2 * side;
  
  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    index_ext = inter->vol_indexR;
    index = inter->vol_indexL;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    index_ext = inter->vol_indexL;
    index = inter->vol_indexR;
  }
  else {
    assert(1==2);
  }


  if (f != NULL){
  
    const unsigned int m = f->model.m;

    //printf("locfa=%d \n",locfa);

    for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {

 
      real flux[m];
      real wL[m];
      for(int iv = 0; iv < m; iv++) {
	wL[iv] = 0;
      }
 
      if (fext != NULL) {  // the right element exists
	real wR[m];
	int ipgR = index_ext[ipgf];
	int ipgL = index[ipgf];
	for(int iv = 0; iv < m; iv++) {
	  int imemR = fext->varindex(fext->deg, fext->raf,fext->model.m, ipgR, iv);
	  wR[iv] = fext->wn[imemR];
	  //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}

	// int_dL F(wL, wR, grad phi_ib)
	real vndsloc[3];
	
	vndsloc[0] = sign * inter->vnds[3 * ipgf + 0];
	vndsloc[1] = sign * inter->vnds[3 * ipgf + 1];
	vndsloc[2] = sign * inter->vnds[3 * ipgf + 2];

	//printf("sign=%d ipgL=%d ipgR=%d vndsloc=%f %f\n",sign,ipgL,ipgR,vndsloc[0],vndsloc[1]);

	f->model.NumFlux(wL, wR, vndsloc, flux);
	//printf("flux=%f %f\n",flux[0],flux[0]);

	// Add flux  to the selected side 
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  f->solver->rhs[imemL] -= flux[iv] * f->dt;
	  printf("imem=%d res=%f\n",imemL,f->solver->rhs[imemL]);
	  
	  /* real* xpg = inter->xpg + 3 * ipgf; */
	  /* real* vnds = inter->vnds + 3 * ipgf; */
	  /* printf("interface flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
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
	  f->solver->rhs[imemL] -= flux[iv] * sign * f->dt;
	  /* printf("boundary flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
	}
      }

    }

  
  }
}

