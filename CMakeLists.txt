set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
cmake_minimum_required(VERSION 2.8)

project(schnaps)

#set(CMAKE_C_COMPILER "tau_cc.sh")
#set(CMAKE_C_FLAGS "-optRevert -optTauSelectFile=select.tau -std=c99 -fopenmp -Ofast")
#set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.8")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g -D_POSIX_C_SOURCE=200112L")

#set(CMAKE_C_FLAGS " -std=c99 -lumfpack -lamd -lblas -lsuitesparseconfig -lcholmod -lcolamd -D_DOUBLE_PRECISION
# -D_SCHNAPS_DIR=${PROJECT_SOURCE_DIR} -O3 -g  ")
#set(CMAKE_C_FLAGS " -std=c99 -D_DOUBLE_PRECISION -D_SCHNAPS_DIR=${PROJECT_SOURCE_DIR} -O3 -g  ")
set(CMAKE_C_FLAGS " -std=gnu99 -D_SCHNAPS_DIR=${PROJECT_SOURCE_DIR} -O3 ")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_DOUBLE_PRECISION")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -g")

# needed for mac
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lsuitesparseconfig")
set(CMAKE_CXX_FLAGS "-lsuitesparseconfig")

MESSAGE(Compil. options: ${CMAKE_C_FLAGS})

OPTION(USE_OPENCL "Activate OpenCL library" ON)

IF(USE_OPENCL) # set with: cmake . -DUSE_OPENCL:BOOL=ON
  MESSAGE(STATUS "Try to find OpenCL Library")
  find_package (OpenCL)
  if(OPENCL_FOUND)
    MESSAGE(STATUS "OpenCL found")
  else()
    MESSAGE(STATUS "OpenCL not found")
  endif(OPENCL_FOUND)
else()
  MESSAGE(STATUS "OpenCL is not activated")
ENDIF(USE_OPENCL)

if(OPENCL_FOUND)
  INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_OPENCL")
else()
  INCLUDE_DIRECTORIES(".")
endif(OPENCL_FOUND)

OPTION(USE_STARPU
  "Compile with the StarPU library"
  ON)

# NB: FindSTARPU.cmake does not properly deal with non-standard
# locations, so it may be necessary to just put the include directory
# in by hand.  For example,
# cmake ../ -DSTARPU_INCDIR=/home/<username>/starpu/include/starpu/1.2/
IF(USE_STARPU)# set with: cmake . -DUSE_STARPU:BOOL=ON
  MESSAGE("Compiling with StarPU Library")
  if(STARPU_INCDIR)
    INCLUDE_DIRECTORIES(${STARPU_INCDIR})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_STARPU")
  else()
    find_package (STARPU REQUIRED COMPONENTS HWLOC MPI)
    if(STARPU_FOUND)
      MESSAGE("StarPU found")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_STARPU")
      INCLUDE_DIRECTORIES(${STARPU_INCLUDE_DIRS})
      LINK_DIRECTORIES(${STARPU_LIBRARY_DIRS})
    else()
      MESSAGE("StarPU not found")
    endif(STARPU_FOUND)
  endif(STARPU_INCDIR)
ENDIF(USE_STARPU)

OPTION(USE_IGRAPH
  "Compile with the igraph library"
  ON)

#INCLUDE_DIRECTORIES("/home/mal/igraph/include/igraph")

# NB: FindSTARPU.cmake does not properly deal with non-standard
# locations, so it may be necessary to just put the include directory
# in by hand.  For example,
# cmake ../ -DIGRAPH_INCDIR=/home/<username>/igraph/include/igraph
IF(USE_IGRAPH)# set with: cmake . -DUSE_IGRAPH:BOOL=ON
  MESSAGE("Compiling with igraph Library")
  if(IGRAPH_INCDIR)
    INCLUDE_DIRECTORIES(${IGRAPH_INCDIR})
  else()
    find_package (igraph)
    if(HAVE_IGRAPH)
      MESSAGE("igraph found")
      INCLUDE_DIRECTORIES(${IGRAPH_INCLUDES})
  else()
    MESSAGE("igraph not found")
  endif(HAVE_IGRAPH)
  endif(IGRAPH_INCDIR)
ENDIF(USE_IGRAPH)

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

#find_package(OpenMP)
#if(OPENMP_FOUND)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

set(SuiteSparse_VERBOSE OFF)
set(SuiteSparse_USE_LAPACK_BLAS ON)
find_package(suitesparse)
if(SuiteSparse_FOUND)
  MESSAGE(STATUS "umfpack libs ok")
  list(APPEND INCLUDE_DIRECTORIES ${SuiteSparse_INCLUDE_DIRS})
  list(APPEND LIBS ${SuiteSparse_LIBRARIES})
  list(APPEND LIBS ${SuiteSparse_LAPACK_BLAS_LIBRARIES})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_UMFPACK")
else()
  MESSAGE(STATUS "umfpack not found: deactivated")
endif()



include_directories(src)

enable_testing()

list(APPEND
  source_files
  src/schnaps.h
  src/simulation.c
  src/simulation.h
  src/simulation_spu.c
  src/simulation_spu.h
  src/linear_solver.h
  src/linear_solver.c
  src/field.c
  src/field.h
  src/interface.h
  src/interface.c
  src/implicit.h
  src/implicit.c
  src/local_implicit.c
  src/geometry.c
  src/geometry.h
  src/global.c
  src/global.h
  src/model.c
  src/model.h
  src/interpolation.c
  src/interpolation.h
  src/macromesh.c
  src/macromesh.h
  src/h20phi.h
  src/maxwell.h
  src/maxwell.c
  src/waterwave2d.h
  src/waterwave2d.c
  src/mhd.h
  src/mhd.c
  src/gyro.h
  src/gyro.c
  src/pic.h
  src/pic.c
  src/skyline.h
  src/skyline.c
  src/skyline_spu.h
  src/skyline_spu.c
  src/csparse.h
  src/csparse.c
  src/csparse_spu.h
  src/csparse_spu.c
  src/dpackfgmres.h
  src/dpackfgmres.c
  src/advanced_linear_solver.h
  src/advanced_linear_solver.c
  src/solvercontinuous.h
  src/solvercontinuous.c
  src/solverpoisson.h
  src/solverpoisson.c
  src/quantities_vp.h
  src/quantities_vp.c
  src/diagnostics_vp.h
  src/diagnostics_vp.c
  src/collision.h
  src/collision.c
  )

if(SuiteSparse_FOUND)
  list(APPEND
    source_files
    src/solverumfpack.h
    src/solverumfpack.c
    )
endif(SuiteSparse_FOUND)


if (OPENCL_FOUND)
  list(APPEND
    source_files
    src/field_cl.c
    src/field_cl.h
    src/clutils.c
    src/clutils.h
    src/clinfo.c
    src/clinfo.h
    )
endif (OPENCL_FOUND)

#schnaps library
add_library (_schnaps ${source_files})

set(LIBS "")

if (OPENCL_FOUND)
  list(APPEND LIBS _schnaps ${OPENCL_LIBRARIES} m)
else()
  list(APPEND LIBS _schnaps  m)
endif (OPENCL_FOUND)


if (STARPU_FOUND)
  list(APPEND LIBS ${STARPU_SHM_LIBRARIES})
  list(APPEND LIBS ${STARPU_LIBRARIES})
  MESSAGE(STATUS "StarPU found at:" ${STARPU_SHM_LIBRARIES})
endif (STARPU_FOUND)

if (HAVE_IGRAPH)
  list(APPEND LIBS  ${IGRAPH_LIBRARIES})
endif (HAVE_IGRAPH)


# nearest neighbour fast search library
find_package(flann)
if(FLANN_FOUND)
  MESSAGE(STATUS "Delicious!")
  list(APPEND LIBS ${FLANN_LIBRARY})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_FLANN")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -D_WITH_FLANN")
endif()

if(SuiteSparse_FOUND)
  list(APPEND LIBS ${SuiteSparse_LIBRARIES})
  list(APPEND LIBS ${SuiteSparse_LAPACK_BLAS_LIBRARIES})
endif()




MESSAGE("detected libs: " ${LIBS})

# main program schnaps
add_executable(
  schnaps
  src/schnaps.c
  )
target_link_libraries(schnaps ${LIBS})


# Find all of the files ending in .c in the test directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "test/")
file(GLOB globtest test/*.c)
MESSAGE(STATUS "globtest")
foreach(test IN LISTS globtest)
  string(REPLACE  ".c"  ""  test  ${test})
  string(REGEX REPLACE  "([^/]+$)"  ""  dir  ${test})
  string(REPLACE  ${dir}  ""  test  ${test})
  #MESSAGE(STATUS ${test})
  list(APPEND testlist ${test})
endforeach(test)

foreach(test IN LISTS testlist)
  MESSAGE(STATUS "Adding test " ${test})
  add_executable(
    ${test}
    test/test.h
    test/${test}.c
    )
  target_link_libraries(${test} ${LIBS})
  add_test(test/${test} test/${test})
endforeach(test)


# Find all of the files ending in .c in the example directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "example/")
file(GLOB globexample example/*.c)
MESSAGE(STATUS "globexample")
foreach(example IN LISTS globexample)
  string(REPLACE  ".c"  ""  example  ${example})
  string(REGEX REPLACE  "([^/]+$)"  ""  dir  ${example})
  string(REPLACE  ${dir}  ""  example  ${example})
  #MESSAGE(STATUS ${example})
  list(APPEND examplelist ${example})
endforeach(example)

foreach(example IN LISTS examplelist)
  MESSAGE(STATUS "Adding example " ${example})
  add_executable(
    ${example}
    test/test.h
    example/${example}.c
    )
  target_link_libraries(${example} ${LIBS})
endforeach(example)
