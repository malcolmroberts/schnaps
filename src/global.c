#include "global.h"

int nplatform_cl = 1;
int ndevice_cl = 0;

bool starpu_is_init = false;
bool starpu_use = false;
bool starpu_c_use = false;
bool starpu_ocl_use = false;

KineticData  schnaps_kinetic_data;

// OpenCL program for StarPU
bool opencl_program_is_init = false;
#ifdef _WITH_STARPU
struct starpu_opencl_program opencl_program;
#endif //_WITH_STARPU

void InitKineticData(KineticData *kd, int nbelemv, int degv){
  kd->nb_elem_v = nbelemv;
  kd->deg_v = degv;
  kd->mv = nbelemv * degv + 1;
  kd->index_max_kin = kd->mv - 1;
  kd->index_rho = kd->mv + 1 ;
  kd->index_phi = kd->mv ;
  kd->index_ex = kd->mv + 2;
  kd->index_ey = kd->mv + 3;
  kd->index_ez = kd->mv + 4;
  kd->index_u = kd->mv + 5;
  kd->index_P = kd->mv + 6;
  kd->index_T = kd->mv + 7;
  kd->index_max = kd->mv + 8;
  kd->vmax = 6;
  kd->dv = 2 * kd->vmax / nbelemv;
  kd->gamma = 3.;
  kd->knud =1.;
}
