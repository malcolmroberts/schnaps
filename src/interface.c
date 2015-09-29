#include "interface.h"
#include <assert.h>
#include <stdio.h>

int VarindexFace(int npg, int m, int ipgf, int iv){

  return ipgf * m + iv;

}

void InitInterface_SPU(Interface* inter){

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
    fd = inter->fL;
    locfa = inter->locfaL;
    npgf = inter->npgL;
    wf = inter->wL;
    wf_handle = inter->wL_handle;
    vol_index = inter->vol_indexL;
    vol_index_handle = inter->vol_indexL_handle;
  }

  if (side == 1) {
    fd = inter->fR;
    locfa = inter->locfaR;
    npgf = inter->npgR;
    wf = inter->wR;
    wf_handle = inter->wR_handle;
    vol_index = inter->vol_indexR;
    vol_index_handle = inter->vol_indexR_handle;
  }
  
  assert(side ==1 || side == 0);

  
  if (fd !=NULL){

    void* arg_buffer;// = malloc(sizeof(int)+sizeof(field));
  size_t arg_buffer_size;
  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, &npgf, sizeof(int),
			   STARPU_VALUE, fd, sizeof(field),
			   0);

    
  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  task->handles[0] = wf_handle;
  task->handles[1] = fd->wn_handle;
  task->handles[2] = vol_index_handle;
  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");


  }
    

}

void ExtractInterface_C(void* buffer[], void* cl_args){

  int npgf;
  field fd;
  
  starpu_codelet_unpack_args(cl_args,
			     &npgf, &fd);
 

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
    for(int iv=0; iv < fd.model.m; iv++){
      int imem = fd.varindex(fd.deg, fd.raf, fd.model.m,
			      ipgv, iv);
      int imemf = VarindexFace(npgf, fd.model.m, ipgf, iv);
      wf[imemf] = wn[imem];
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

void InterfaceExplicitFlux_SPU(Interface* inter, int side)
{


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

