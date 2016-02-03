#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int TestCodelet_ZeroBuffer_SPU(void){
  bool test = true;


  // Kernel compilation options
  sprintf(cl_buildoptions, "%s", "");
  char buf[1000];
#ifdef _DOUBLE_PRECISION
  sprintf(buf, "-D real=double");
#else
  sprintf(buf, "-D real=float");
#endif
  strcat(cl_buildoptions, buf);


  // Data buffer
  const int size = 10000000;
  schnaps_real* buffer = calloc(size, sizeof(schnaps_real));

  putenv("STARPU_NCUDA=0");
  // Init StarPU
  test &= (starpu_init(NULL) != -ENODEV);

  // Activate all the available codelets
  starpu_c_use = true;
  starpu_ocl_use = true;

  // Interesting environment variables
  printf("Environment variables...\n");
  printf("STARPU_NCPU                : %s\n", getenv("STARPU_NCPU"));
  printf("STARPU_NOPENCL             : %s\n", getenv("STARPU_NOPENCL"));
  printf("STARPU_NCUDA               : %s\n", getenv("STARPU_NCUDA"));
  printf("STARPU_NMIC                : %s\n", getenv("STARPU_NMIC"));
  printf("STARPU_OPENCL_ON_CPUS      : %s\n", getenv("STARPU_OPENCL_ON_CPUS"));
  printf("STARPU_OPENCL_ONLY_ON_CPUS : %s\n", getenv("STARPU_OPENCL_ONLY_ON_CPUS"));

  // Get workers
  printf("Available workers...\n");
  const int nb_workers = starpu_worker_get_count();
  printf("Number of workers          : %d\n", nb_workers);
  printf("CPU workers                : %d\n", starpu_cpu_worker_get_count());
  const int nb_ocl = starpu_opencl_worker_get_count();
  printf("OCL workers                : %d\n", nb_ocl);
  printf("CUDA workers               : %d\n", starpu_cuda_worker_get_count());
  printf("MIC workers                : %d\n", starpu_mic_worker_get_count());
  
  // Loop over workers to submit tasks
  for (int wid = 0; wid < nb_workers; ++wid) {
    // Create context with a single worker
    unsigned int ctxid = starpu_sched_ctx_create(&wid, 1, "ctx", STARPU_SCHED_CTX_POLICY_NAME, "prio", NULL);
    starpu_sched_ctx_set_context(&ctxid);
  
    // Create data handle (init and register)
    for (int i = 0; i < size; ++i) buffer[i] = i;
    starpu_data_handle_t handle;
    starpu_vector_data_register(&handle, 0, (uintptr_t) buffer, size, sizeof(schnaps_real));

    // Task
    ZeroBuffer_SPU(handle);

    // Check output
    starpu_task_wait_for_all();
    starpu_data_prefetch_on_node(handle, 0, 0);
    starpu_data_unregister(handle);
    for (int i = 0; i < size; ++i) test &= (abs(buffer[i]) < _VERY_SMALL);
    
    // Delete context
    starpu_sched_ctx_delete(ctxid);
  }

  // Delete opencl program if it has been created
  if (nb_ocl > 0) {
    int ret = starpu_opencl_unload_opencl(&opencl_program);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_unload_opencl");
  }

  starpu_shutdown();


  return test;
}

int main(void) {
  // Unit tests
  int resu = TestCodelet_ZeroBuffer_SPU();
  if (resu) printf("StarPU ZeroBuffer Codelet test OK !\n");
  else printf("StarPU ZeroBuffer Codelet test failed !\n");
  return !resu;
}
