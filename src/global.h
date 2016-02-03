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
#endif //_WITH_STARPU


#endif // #ifndef _GLOBAL_H
