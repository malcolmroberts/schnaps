#include "clinfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

void InitCLInfo(CLInfo* cli,int platform_id,int device_id){
  cl_int status;  // for checking OpenCL errors

  // numbers of platforms
  status = clGetPlatformIDs(0, NULL, &(cli->nbplatforms));
  assert(status == CL_SUCCESS);

  assert(cli->nbplatforms > 0);

  // platform array construction
  cl_platform_id* platforms = malloc(sizeof(cl_platform_id[cli->nbplatforms]));
  assert(platforms);

  status = clGetPlatformIDs(cli->nbplatforms, platforms, NULL);
  assert(status == CL_SUCCESS);

  printf("Available platforms:\n");

  char pbuf[2000];
  for(int i = 0; i < cli->nbplatforms; ++i) {
    printf("\nPlatform %d:\n", i);
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_NAME,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    assert(status == CL_SUCCESS);
    printf("%s\n", pbuf);

    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VENDOR,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    assert(status == CL_SUCCESS);
    printf("%s\n",pbuf);

    //  opencl version
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VERSION,
			       sizeof(cli->platformname),
			       cli->platformname,
			       NULL);
    assert(status == CL_SUCCESS);
    printf("%s\n",cli->platformname);
  }

  cli->platformid=platform_id;
  assert(cli->platformid < cli->nbplatforms);

  // devices count
  status = clGetDeviceIDs(platforms[platform_id],
			  CL_DEVICE_TYPE_ALL,
			  0,
			  NULL,
			  &(cli->nbdevices));
  assert(status == CL_SUCCESS);
  assert(cli->nbdevices > 0);

  // devices array construction
  cli->device = malloc(sizeof(cl_device_id) * cli->nbdevices);
  assert(cli->device);


  printf("\nWe choose device %d/%d ",device_id,cli->nbdevices-1);
  assert(device_id < cli->nbdevices);
  printf("of platform %d/%d\n",platform_id,cli->nbplatforms-1);
  status = clGetDeviceIDs(platforms[platform_id],
			  CL_DEVICE_TYPE_ALL,
			  cli->nbdevices,
			  cli->device,
			  NULL);
  assert(status == CL_SUCCESS);

  // device name
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_NAME,
			   sizeof(cli->devicename),
			   cli->devicename,
			   NULL);
  assert(status == CL_SUCCESS);
  printf("%s\n",cli->devicename);


  // device memory
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_GLOBAL_MEM_SIZE,
			   sizeof(cl_ulong),
			   &(cli->devicememsize),
			   NULL);
  assert(status == CL_SUCCESS);
  printf("Global memory: %f MB\n",cli->devicememsize/1024./1024.);

  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_MEM_ALLOC_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxmembuffer),
			   NULL);
  assert(status == CL_SUCCESS);
  printf("Max buffer size: %f MB\n",cli->maxmembuffer/1024./1024.);

  // compute unit size cache
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_LOCAL_MEM_SIZE,
			   sizeof(cl_ulong),
                           &(cli->cachesize),
			   NULL);
  assert(status == CL_SUCCESS);
  printf("Cache size: %f KB\n",cli->cachesize/1024.);

  // get maxconstmem
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
			   sizeof(cl_ulong),
                           &(cli->maxconstmem),
			   NULL);
  assert(status == CL_SUCCESS);
  printf("Const mem: %f KB\n",cli->maxconstmem/1024.);

  // get maxconst args
  int maxcstargs;
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_CONSTANT_ARGS,
			   sizeof(cl_ulong),
                           &maxcstargs,
			   NULL);
  assert(status == CL_SUCCESS);
  printf("Max Const args: %d \n",maxcstargs);

  // nb of compute units
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_COMPUTE_UNITS,
			   sizeof(cl_uint),
			   (void*)&(cli->nbcomputeunits),
			   NULL);
  assert(status == CL_SUCCESS);
  printf("Nb of compute units: %d\n",cli->nbcomputeunits);

  // max workgroup size
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_WORK_GROUP_SIZE,
			   sizeof(size_t),
			   (void*)&(cli->maxworkgroupsize),
			   NULL);
  assert(status == CL_SUCCESS);
  printf("Max workgroup size: %zu\n",cli->maxworkgroupsize);

  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_EXTENSIONS,
			   sizeof(cli->clextensions),
			   cli->clextensions,
			   NULL);
  assert(status == CL_SUCCESS);
  printf("OpenCL extensions:\n%s\n",cli->clextensions);

  // first opencl context
  cli->context = clCreateContext(NULL, // no context properties
				 1,         // only one device in the list
				 &cli->device[device_id], // device list
				 NULL, // callback function
				 NULL, // function arguments
				 &status);
  assert(status == CL_SUCCESS);

  // command queue
  cli->commandqueue = clCreateCommandQueue(cli->context,
					   cli->device[device_id],
					   CL_QUEUE_PROFILING_ENABLE,
					   &status);
  assert(status == CL_SUCCESS);

  printf("Init OK\n\n");
}

void PrintCLInfo(CLInfo* cli){
  cl_int status;  // for checking OpenCL errors
  char pbuf[2000];

  printf("%s\n",cli->platformname);
  printf("%s\n",cli->devicename);

  // device memory
  printf("Global memory: %f MB\n",cli->devicememsize/1024./1024.);
  printf("Max buffer size: %f MB\n",cli->maxmembuffer/1024./1024.);
  printf("Cache size: %f KB\n",cli->cachesize/1024.);
  printf("Nb of compute units: %d\n",cli->nbcomputeunits);
  printf("Max workgroup size: %zu\n",cli->maxworkgroupsize);
  printf("OpenCL extensions:\n%s\n",cli->clextensions);
}

void BuildKernels(CLInfo* cli,char* strprog){
  cl_int err;

  // kernels creation
  cli->program =
    clCreateProgramWithSource(cli->context,
                              1,
			      (const char **) &strprog,
                              NULL,
                              &err);
  //printf("%s\n",strprog);
  //assert(1==2);
  if(!(cli->program)) {
    printf("Failed to create program.\n");
  }

  // compilation
  err = clBuildProgram(cli->program,
		       0,               // one device
		       NULL,
		       NULL, NULL, NULL);

  /* cl_int clBuildProgram(	cl_program program, */
  /* 				cl_uint num_devices, */
  /* 				const cl_device_id *device_list, */
  /* 				const char *options, */
  /* 				void(CL_CALLBACK *pfn_notify)(cl_program program, void *user_data), */
  /* 				void *user_data) */

  // if not successfull: display the errors
  if(err != CL_SUCCESS) {
    size_t len;
    char buffer[1024 * 1024];
    printf("Compilation output:\n");
    clGetProgramBuildInfo(cli->program,
			  cli->device[cli->deviceid],
			  CL_PROGRAM_BUILD_LOG, sizeof(buffer),
			  buffer,
			  &len);
    printf("len=%d buf=%zd\n",(int) len,sizeof(buffer));
    printf("%s\n",buffer);
  }
  assert( err == CL_SUCCESS);
}

void ReadFile(char filename[],char** s){
  FILE* f = fopen( filename , "r" );
  assert(f);

  fseek( f , 0L , SEEK_END);
  int size = ftell(f);
  rewind(f);

  /* allocate memory for entire content */
  *s = calloc(size + 2, sizeof(char));
  assert(*s);
  /* copy the file into the buffer */
  int status = fread(*s , sizeof(char), size , f);
  // assert(status==0);
  //assert(s[size]=='\0');
  fclose(f);
}


//! \brief scan all *.h and *.c in order to find the code
//! to be shared with opencl
//! such code is enclosed between #pragma start_opencl and 
//! #pragma end_opencl
void GetOpenCLCode(char *model){

  int status;
  char getModel[50];
  strcpy(getModel, "sh get_model.sh ");
  strcat(getModel, model);

  printf("ligne de commande : %s\n",getModel);

  status = system(getModel);

  status=system("sh get_opencl_code.sh");
  assert(!status);
 
}
  
