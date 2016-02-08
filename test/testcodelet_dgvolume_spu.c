#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

bool submit_task(Simulation* simu, schnaps_real* buffer_cmp) {
  bool test = true;

  // Indicators for function comparision
  static bool first_function = true;
  static bool second_function = false;

  // Data buffer
  const int size = simu->wsize / simu->macromesh.nbelems;
  schnaps_real* buffer_in = calloc(size, sizeof(schnaps_real));
  schnaps_real* buffer_out = calloc(size, sizeof(schnaps_real));

  // Loop over the fields
  for (int iel = 0; iel < simu->macromesh.nbelems; ++iel) {
    field* f = simu->fd + iel;

    // Create data handle (init and register)
    for (int i = 0; i < size; ++i) buffer_in[i] = i;
    starpu_vector_data_register(&f->wn_handle, 0, (uintptr_t) buffer_in, size, sizeof(schnaps_real));
    starpu_vector_data_register(&f->res_handle, 0, (uintptr_t) buffer_out, size, sizeof(schnaps_real));

    // Task
    DGVolume_SPU(f);

    // Check output
    starpu_task_wait_for_all();
    starpu_data_prefetch_on_node(f->wn_handle, 0, 0);
    starpu_data_unregister(f->wn_handle);
    starpu_data_prefetch_on_node(f->res_handle, 0, 0);
    starpu_data_unregister(f->res_handle);

    for (int i = 0; i < size; ++i)
      assert(abs(buffer_in[i] - i) < _VERY_SMALL);

    // Store the first call, otherwise compare
    if (first_function) {
      for (int i = 0; i < size; ++i)
        buffer_cmp[iel * size + i] = buffer_out[i];
    } else {
      for (int i = 0; i < size; ++i) {
        //if (!(abs(buffer_cmp[iel * size + i] - buffer_out[i]) < _VERY_SMALL))
        //printf("field: %d\twid: %d\tocl: %f\tcpu: %f\n",
        //       iel, i, buffer_cmp[iel * size + i], buffer_out[i]);
        test &= (abs(buffer_cmp[iel * size + i] - buffer_out[i]) < _VERY_SMALL);
      }
    }
  }

  if (first_function) {
    first_function = false;
    second_function = true;
    printf(" result set to comparision reference.\n");
  } else {
    if (test) printf(" OK\n");
    else printf(" KO !\n");

    if (second_function && !test)
      printf("WARNING: Error may come from first function.\n");
    if (second_function)
      second_function = false;
  }

  free(buffer_in);
  free(buffer_out);

  return test;
}


int TestCodelet_DGVolume_SPU(void){
  bool test = true;

  /* int deg[]={3, 3, 3}; */
  /* int raf[]={3, 3, 3}; */
  int deg[]={1, 2, 2};
  int raf[]={2, 2, 2};

  MacroMesh mesh;
  //ReadMacroMesh(&mesh,"../test/testdisque.msh");
  ReadMacroMesh(&mesh,"../test/testcube2.msh");
  BuildConnectivity(&mesh);
  CheckMacroMesh(&mesh, deg, raf);

  Model model;
  model.m = 6;
  model.NumFlux = Maxwell3DNumFlux_upwind;
  model.BoundaryFlux = Maxwell3DBoundaryFlux_upwind;
  model.InitData = Maxwell3DInitData;
  model.ImposedData = Maxwell3DImposedData;
  model.Source = NULL;

  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);


  // Kernel compilation options
  sprintf(cl_buildoptions, "%s", "");
  char buf[1000];
#ifdef _DOUBLE_PRECISION
  sprintf(buf, "-D schnaps_real=double");
#else
  sprintf(buf, "-D schnaps_real=float");
#endif
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D NUMFLUX=%s", "Maxwell3DNumFlux_upwind");
  strcat(cl_buildoptions, buf);


  // Init StarPU
  struct starpu_conf conf;
  test &= (starpu_conf_init(&conf) == 0);
  assert(test);
  test &= (starpu_init(&conf) != -ENODEV);
  assert(test);

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

  // Workers
  const int nb_workers = starpu_worker_get_count();
  const int nb_ocl = starpu_opencl_worker_get_count();
  printf("Available workers...\n");
  printf("Number of workers          : %d\n", nb_workers);
  printf("CPU workers                : %d\n", starpu_cpu_worker_get_count());
  printf("OPENCL workers             : %d\n", nb_ocl);
  printf("CUDA workers               : %d\n", starpu_cuda_worker_get_count());
  printf("MIC workers                : %d\n", starpu_mic_worker_get_count());

  // Regenerate opencl code
  if (nb_ocl > 0) GetOpenCLCode();

  // Codelet
  struct starpu_codelet* codelet = DGVolume_codelet();
  struct starpu_codelet codelet_backup = *codelet;

  schnaps_real* buffer_cmp = calloc(simu.wsize, sizeof(schnaps_real));

  // Loop over workers to submit tasks
  for (int wid = 0; wid < nb_workers; ++wid) {
    // Create context with a single worker
    unsigned int ctxid = starpu_sched_ctx_create(
        &wid, 1, "ctx",
        // TRICK Default scheduler:
        // avoids seg fault in _starpu_push_task_on_specific_worker
        STARPU_SCHED_CTX_POLICY_NAME, conf.sched_policy_name,
        NULL);
    starpu_sched_ctx_set_context(&ctxid);
    starpu_sched_ctx_display_workers(ctxid, stdout);

    // Empty codelet for function selection
    for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
      codelet->cpu_funcs[i] = NULL;
      codelet->opencl_funcs[i] = NULL;
      codelet->cuda_funcs[i] = NULL;
      codelet->mic_funcs[i] = NULL;
    }

    // TRICK Need at least one C function (at least when OpenCL worker)
    codelet->cpu_funcs[0] = codelet_backup.cpu_funcs[0];

    // Loop over the codelet implementations
    switch (starpu_worker_get_type(wid)) {
      case STARPU_CPU_WORKER:
        if (codelet_backup.cpu_funcs[0] == NULL) {
          printf("No C codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.cpu_funcs[i] != NULL) {
              codelet->cpu_funcs[0] = codelet_backup.cpu_funcs[i];
              printf("Submit C codelet %d...", i);
              test &= submit_task(&simu, buffer_cmp);
            }
          }
        }
        break;

      case STARPU_OPENCL_WORKER:
        if (codelet_backup.opencl_funcs[0] == NULL) {
          printf("No OpenCL codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.opencl_funcs[i] != NULL) {
              codelet->opencl_funcs[0] = codelet_backup.opencl_funcs[i];
              printf("Submit OpenCL codelet %d...", i);
              test &= submit_task(&simu, buffer_cmp);
            }
          }
        }
        break;

      case STARPU_CUDA_WORKER:
        if (codelet_backup.cuda_funcs[0] == NULL) {
          printf("No CUDA codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.cuda_funcs[i] != NULL) {
              codelet->cuda_funcs[0] = codelet_backup.cuda_funcs[i];
              printf("Submit CUDA codelet %d...", i);
              test &= submit_task(&simu, buffer_cmp);
            }
          }
        }
        break;

      case STARPU_MIC_WORKER:
        if (codelet_backup.mic_funcs[0] == NULL) {
          printf("No MIC codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.mic_funcs[i] != NULL) {
              codelet->mic_funcs[0] = codelet_backup.mic_funcs[i];
              printf("Submit MIC codelet %d...", i);
              test &= submit_task(&simu, buffer_cmp);
            }
          }
        }
        break;

      default:
        printf("Untreated worker type.\n");
    }

    // Delete context
    starpu_sched_ctx_delete(ctxid);
  }

  *codelet = codelet_backup;

  // Delete opencl program if it has been created
  if (nb_ocl > 0) {
    int ret = starpu_opencl_unload_opencl(&opencl_program);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_unload_opencl");
  }

  starpu_shutdown();

  free(buffer_cmp);

  FreeMacroMesh(&mesh);


  return test;
}

int main(void) {
  // Unit tests
  int resu = TestCodelet_DGVolume_SPU();
  if (resu) printf("StarPU DGVolume Codelet test OK !\n");
  else printf("StarPU DGVolume Codelet test failed !\n");
  return !resu;
}
