#ifndef _GLOBAL_H
#define _GLOBAL_H
//! brief global variables and defs

#include <stdbool.h>

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
#define real float
#else
#define real double
#endif

#ifndef _DOUBLE_PRECISION
#define  _VERY_SMALL (1e-8)
#define _SMALL (1e-5)
#else
#define _VERY_SMALL (1e-16)
#define _SMALL (1e-10)
#endif


bool starpu_is_init;
bool starpu_use;


#endif // #ifndef _GLOBAL_H
