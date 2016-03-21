#ifndef _GLOBAL_H
#define _GLOBAL_H
//! brief global variables and defs

#include <stdbool.h>
#ifdef _WITH_STARPU
#include <starpu.h>
#endif //_WITH_STARPU

#ifndef _OPENMP
// activate pthread if openmp is not here
//#define _WITH_PTHREAD
#endif

#ifdef _WITH_OPENCL
extern int nplatform_cl;
extern int ndevice_cl;
char numflux_cl_name[1024]; // FIXME: move to field struct.
char cl_buildoptions[1024];
#endif //_WITH_OPENCL

#define __constant
#define __local
#define __private

#ifndef _DOUBLE_PRECISION
#define schnaps_real float
#else
#define schnaps_real double
#endif

#ifndef _DOUBLE_PRECISION
#define  _VERY_SMALL (1e-8)
#define _SMALL (1e-5)
#else
#define _VERY_SMALL (1e-16)
#define _SMALL (1e-10)
#endif


typedef struct KineticData{
  int time_order;
  int nb_elem_v;
  int deg_v;
  int mv;
  int index_max_kin;
  int index_rho;
  int index_phi;
  int index_ex;
  int index_ey;
  int index_ez;
  int index_u;
  int index_P;
  int index_T;
  int index_max;
  schnaps_real vmax;
  schnaps_real dv;
  schnaps_real gamma;
  schnaps_real knud;
  bool solve_quasineutrality;
  bool substract_mean_charge;
  // quasi neutrality damping (term in front of (phi -phibar) in
  // QN or damped Laplace equation)
  schnaps_real qn_damping;
} KineticData;


typedef struct LatticeData{
  int q;
  int d;
  double ** q_tab;
  double * w_tab;
  
  int temp_const;
  double c;

  int index_max_q;
  int index_rho;
  int index_ux;
  int index_uy;
  int index_uz;
  int index_temp;
  int index_p;
  int index_max;
  
} LatticeData;


extern KineticData schnaps_kinetic_data;
extern LatticeData schnaps_lattice_data;

void InitKineticData(KineticData *kd, int nbelemv, int degv);
void InitLatticeData(LatticeData *kd, int dim, int Q,int temp,double sound);

extern bool starpu_is_init;
extern bool starpu_use;
// Use of c kernels when use of starpu
extern bool starpu_c_use;
// Use of opencl kernels when use of starpu
extern bool starpu_ocl_use;

// OpenCL program for StarPU
extern bool opencl_program_is_init;
#ifdef _WITH_STARPU
extern struct starpu_opencl_program opencl_program;

#define LOAD_OPENCL_PROGRAM_SPU()               \
  if (!opencl_program_is_init) {                \
    opencl_program_is_init = true;              \
    printf("load OpenCL program...\n");         \
    STARPU_CHECK_RETURN_VALUE(                  \
        starpu_opencl_load_opencl_from_file(    \
            "./schnaps.cl",                     \
            &opencl_program,                    \
            cl_buildoptions),                   \
        "starpu_opencl_load_opencl_from_file"); \
  }

#else

#define LOAD_OPENCL_PROGRAM_SPU()               \
  assert(opencl_program_is_init);

#endif //_WITH_STARPU



#endif // #ifndef _GLOBAL_H
