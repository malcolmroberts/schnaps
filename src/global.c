#include "global.h"

int nplatform_cl = 0;
int ndevice_cl = 1;

bool starpu_is_init = false;
bool starpu_use = false;
bool starpu_c_use = false;
bool starpu_ocl_use = false;

// OpenCL program for StarPU
bool opencl_program_is_init = false;
#ifdef _WITH_STARPU
struct starpu_opencl_program opencl_program;
#endif //_WITH_STARPU
