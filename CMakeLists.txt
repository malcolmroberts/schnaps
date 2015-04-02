
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
cmake_minimum_required(VERSION 2.8)

project(schnaps)

#set(CMAKE_C_COMPILER "tau_cc.sh")
#set(CMAKE_C_FLAGS "-optRevert -optTauSelectFile=select.tau -std=c99 -fopenmp -Ofast")
#set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.8")

#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -g")
set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3")
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()



include_directories(src)

enable_testing()

list(APPEND
source_files
src/field.h
src/schnaps.h
src/geometry.h
src/global.h
src/h20phi.h
src/interpolation.h
src/macromesh.h
src/model.h
src/collision.h
src/gyro.h
src/skyline.h
src/solverpoisson.h
src/quantities_vp.h
src/diagnostics_vp.h
src/field.c
src/geometry.c
src/interpolation.c
src/macromesh.c
src/model.c
src/collision.c
src/gyro.c
src/skyline.c
src/solverpoisson.c
src/quantities_vp.c
src/diagnostics_vp.c
)

#schnaps library
add_library (_schnaps ${source_files})

# main program schnaps
add_executable(
  schnaps
  src/schnaps.c
)
target_link_libraries(schnaps _schnaps m pthread)

# list of tests
add_executable(
  testmacromesh
  test/test.h
  test/test_mesh/testmacromesh.c
)
target_link_libraries(testmacromesh _schnaps m pthread)
add_test(testmacromesh testmacromesh)

add_executable(
  testgeometry
  test/test.h
  test/test_mesh/testgeometry.c
)
target_link_libraries(testgeometry _schnaps m pthread)
add_test(testgeometry testgeometry)

add_executable(
  testinterpolation
  test/test.h
  test/testinterpolation.c
)
target_link_libraries(testinterpolation _schnaps m pthread)
add_test(testinterpolation testinterpolation)

add_executable(
  testmodel
  test/test.h
  test/testmodel.c
)
target_link_libraries(testmodel _schnaps m pthread)

add_test(testfield testfield)
add_executable(
  testfield
  test/test.h
  test/test_field/testfield.c
)
target_link_libraries(testfield _schnaps m pthread)
add_test(testfield testfield)

add_test(testfielddg testfielddg)
add_executable(
  testfielddg
  test/test.h
  test/test_field/testfielddg.c
)
target_link_libraries(testfielddg _schnaps m pthread)
add_test(testfielddg testfielddg)

add_test(testfieldsource testfieldsource)
add_executable(
  testfieldsource
  test/test.h
  test/test_field/testfieldsource.c
)
target_link_libraries(testfieldsource _schnaps m pthread)
add_test(testfieldsource testfieldsource)
 
add_test(testfieldrk2 testfieldrk2)
add_executable(
  testfieldrk2
  test/test.h
  test/test_field/testfieldrk2.c
)
target_link_libraries(testfieldrk2 _schnaps m pthread)
add_test(testfieldrk2 testfieldrk2)

add_test(test2dmeshdetection test2dmeshdetection)
add_executable(
  test2dmeshdetection
  test/test.h
  test/test_mesh/test2dmeshdetection.c
)
target_link_libraries(test2dmeshdetection _schnaps m pthread)
add_test(test2dmeshdetection test2dmeshdetection)

add_test(testfieldrk2-2d testfieldrk2-2d)
add_executable(
  testfieldrk2-2d
  test/test.h
  test/test_field/testfieldrk2-2d.c
)
target_link_libraries(testfieldrk2-2d _schnaps m pthread)
add_test(testfieldrk2-2d testfieldrk2-2d)


add_test(testsubcelldgvol testsubcelldgvol)
add_executable(
  testsubcelldgvol
  test/test.h
  test/test_mesh/testsubcelldgvol.c
)
target_link_libraries(testsubcelldgvol _schnaps m pthread)
add_test(testsubcelldgvol testsubcelldgvol)

add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)
add_executable(
  testfieldrk2-2d-subcell
  test/test.h
  test/test_field/testfieldrk2-2d-subcell.c
)
target_link_libraries(testfieldrk2-2d-subcell _schnaps m pthread)
add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)

add_test(testmeq2 testmeq2)
add_executable(
  testmeq2
  test/test.h
  test/testmeq2.c
)
target_link_libraries(testmeq2 _schnaps m pthread)
add_test(testmeq2 testmeq2)

add_test(testgyro testgyro)
add_executable(
  testgyro
  test/test.h
  test/testgyro.c
)
target_link_libraries(testgyro _schnaps m pthread)
add_test(testgyro testgyro)

add_test(testperiodic testperiodic)
add_executable(
  testperiodic
  test/test.h
  test/test_vp/testperiodic.c
)
target_link_libraries(testperiodic _schnaps m pthread)
add_test(testperiodic testperiodic)

add_test(testskyline testskyline)
add_executable(
  testskyline
  test/test.h
  test/test_vp/testskyline.c
  )
target_link_libraries(testskyline _schnaps m pthread)
add_test(testskyline testskyline)
 
add_test(testpoisson testpoisson)
add_executable(
  testpoisson
  test/test.h
  test/test_vp/testpoisson.c
)
target_link_libraries(testpoisson _schnaps m pthread)
add_test(testpoisson testpoisson)

add_test(test_transportvp test_transportvp)
add_executable(
  test_transportvp
  test/test.h
  test/test_vp/test_transportvp.c
)
target_link_libraries(test_transportvp _schnaps m pthread)
add_test(test_transportvp test_transportvp)

