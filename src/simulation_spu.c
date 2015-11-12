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



/* void DtFields_SPU(Simulation *simu, */
/* 		  starpu_data_handle_t* w_handle, */
/* 		  starpu_data_handle_t* dtw_handle) */
/* { */

/*   /\*   if(simu->pre_dtfields != NULL) { *\/ */
/*   /\*   simu->pre_dtfields(simu, w); *\/ */
/*   /\*   //assert(1==2); *\/ */
/*   /\* } *\/ */

/* #ifdef _OPENMP */
/* #pragma omp parallel */
/* #endif */

/*   //real *w = simu->fd[0].wn; */
/*   //real *dtw = simu->fd[0].dtwn; */

/*   int fsize =  simu->wsize / simu->macromesh.nbelems; */

/*   for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) ZeroBuffer_SPU(dtw_handle[ie]); */

/*   for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) { */
/*     simu->fd[ie].tnow = simu->tnow; */
/*   }  */

/*   for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){ */
/*     Interface* inter = simu->interface + ifa; */
/*     // left = 0  right = 1 */
/*     ExtractInterface_SPU(inter, 0); */
/*     ExtractInterface_SPU(inter, 1); */
/*     if (inter->fR != NULL) { */
/*       InterfaceExplicitFlux_SPU(inter, 0); */
/*       InterfaceExplicitFlux_SPU(inter, 1); */
/*     } */
/*     else{ */
/*       InterfaceBoundaryFlux_SPU(inter); */
/*     } */
/*   } */

  
/*   for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) { */
/*     DGSubCellInterface(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */
/*     DGVolume(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */
/*     DGSource(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */
/*     DGMass(simu->fd + ie, w + ie * fsize, dtw + ie * fsize); */

/*   } */

/*   if(simu->post_dtfields != NULL) { */
/*     simu->post_dtfields(simu, w); */
/*     //assert(1==2); */
/*   } */



/* } */
