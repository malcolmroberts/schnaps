set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
cmake_minimum_required(VERSION 2.8)

project(schnaps)

#set(CMAKE_C_COMPILER "tau_cc.sh")
#set(CMAKE_C_FLAGS "-optRevert -optTauSelectFile=select.tau -std=c99 -fopenmp -Ofast")
#set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.8")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g -D_POSIX_C_SOURCE=200112L")
#set(CMAKE_C_FLAGS "-std=gnu99 -D_REENTRANT -O3 -g")

#set(CMAKE_C_FLAGS " -std=c99 -D_REENTRANT -O3 ")
set(CMAKE_C_FLAGS " -std=c99  -D_DOUBLE_PRECISION")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_DOUBLE_PRECISION")

OPTION(USE_OPENCL
  "Compile with the OpenCL library"
  OFF)

IF(USE_OPENCL)# set with: cmake . -DUSE_OPENCL:BOOL=ON
  MESSAGE(STATUS "Compiling with OpenCL Library")
  find_package (OpenCL)
  if(OPENCL_FOUND)
    MESSAGE(STATUS "OpenCL found")
  else()
    MESSAGE(STATUS "OpenCL not found")
  endif(OPENCL_FOUND)
ENDIF(USE_OPENCL)

if(OPENCL_FOUND)
  INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_OPENCL")
else()
  INCLUDE_DIRECTORIES(".")
endif(OPENCL_FOUND)

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET(BUILD_PARALUTION      OFF CACHE BOOL   "Build paralution library")

IF(BUILD_PARALUTION)

  ADD_DEFINITIONS(-DPARALUTION)
  INCLUDE(ExternalProject)
  SET(PARALUTION_VERSION 1.0.0 CACHE STRING "Paralution version number")
  EXTERNALPROJECT_ADD( paralution-${PARALUTION_VERSION}
     URL  http://www.paralution.com/downloads/paralution-${PARALUTION_VERSION}.tar.gz
     CMAKE_ARGS -DSUPPORT_CUDA=OFF -DSUPPORT_OMP=OFF
   -DBUILD_EXAMPLES=OFF -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -O3
     INSTALL_COMMAND ""
  )

  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/paralution-${PARALUTION_VERSION}-prefix/src/paralution-${PARALUTION_VERSION}/src/)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/paralution-${PARALUTION_VERSION}-prefix/src/paralution-${PARALUTION_VERSION}-build/lib)
  ADD_LIBRARY(paralution_c   src/paralution_c.cpp )
  ADD_DEPENDENCIES(paralution_c  paralution-${PARALUTION_VERSION})
  TARGET_LINK_LIBRARIES(paralution_c paralution)

ENDIF()


include_directories(src)

enable_testing()

list(APPEND
  source_files
  src/schnaps.h
  src/simulation.c
  src/simulation.h
  src/field.c
  src/field.h
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
  src/dpackfgmres.h
  src/dpackfgmres.c
  src/linear_solver.h
  src/linear_solver.c
  src/solverpoisson.h
  src/solverpoisson.c
  src/quantities_vp.h
  src/quantities_vp.c
  src/diagnostics_vp.h
  src/diagnostics_vp.c
  src/collision.h
  src/collision.c
  )


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

IF(BUILD_PARALUTION)
  target_link_libraries(_schnaps paralution_c paralution)
ENDIF(BUILD_PARALUTION)

set(LIBS "")

if (OPENCL_FOUND)
  list(APPEND LIBS _schnaps pthread ${OPENCL_LIBRARIES} m)
else()
  list(APPEND LIBS _schnaps pthread m)
endif (OPENCL_FOUND)

# nearest neighbour fast search library
#find_package(flann)
if(FLANN_FOUND)
  MESSAGE(STATUS "Delicious!")
  list(APPEND LIBS ${FLANN_LIBRARY})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_FLANN")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -D_WITH_FLANN")
endif()


# list of tests
add_executable(
  testmacromesh
  test/test.h
  test/testmacromesh.c
  )
target_link_libraries(testmacromesh ${LIBS})
add_test(testmacromesh testmacromesh)

add_executable(
  testgeometry
  test/test.h
  test/testgeometry.c
  )
target_link_libraries(testgeometry ${LIBS})
add_test(testgeometry testgeometry)

add_executable(
  testinterpolation
  test/test.h
  test/testinterpolation.c
  )
target_link_libraries(testinterpolation ${LIBS})
add_test(testinterpolation testinterpolation)

add_executable(
  testmodel
  test/test.h
  test/testmodel.c
  )
target_link_libraries(testmodel ${LIBS})

add_test(testfield testfield)
add_executable(
  testfield
  test/test.h
  test/testfield.c
  )
target_link_libraries(testfield ${LIBS})
add_test(testfield testfield)

add_test(testfielddg testfielddg)
add_executable(
  testfielddg
  test/test.h
  test/testfielddg.c
  )
target_link_libraries(testfielddg ${LIBS})
add_test(testfielddg testfielddg)

#add_test(testfieldrk2 testfieldrk2)
#add_executable(
#  testfieldrk2
#  test/test.h
#  test/testfieldrk2.c
#  )
#target_link_libraries(testfieldrk2 ${LIBS})
#add_test(testfieldrk2 testfieldrk2)


