cmake_minimum_required(VERSION 3.19)
project(COCOGoldenTest
        DESCRIPTION "COCO golden master tests"
        LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Disable floating-point contraction (FMA), see src/coco_fp_contract.h
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-ffp-contract=off COCO_HAVE_FP_CONTRACT_OFF)
if(COCO_HAVE_FP_CONTRACT_OFF)
  add_compile_options(-ffp-contract=off)
endif()

include(CTest)

find_library(LIBM m)

## The library is built separately from the driver, which is written in C, so
## that the driver only sees the public C interface.
add_library(coco STATIC coco.cpp coco.h)
target_include_directories(coco PUBLIC .)
if(LIBM)
  target_link_libraries(coco PUBLIC ${LIBM})
endif()
## With optimization enabled, clang replaces pow(10, x) by exp10(x) and
## pow(2, x) by exp2(x), which can change results in the last bit. Always
## optimize like the Python package does (meson debugoptimized) so that the
## results do not depend on the build type.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(coco PRIVATE -O2)
endif()

add_executable(golden golden.c)
target_link_libraries(golden PRIVATE coco)

## Results depend on the platform's libm, so references are kept per platform.
set(GOLDEN_PLATFORM "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
set(GOLDEN_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/reference/${GOLDEN_PLATFORM}")

set(GOLDEN_SUITES
  toy
  bbob
  bbob-boxed
  bbob-noisy
  bbob-largescale
  bbob-mixint
  bbob-biobj
  bbob-biobj-ext
  bbob-biobj-mixint
  bbob-constrained
  bbob-constrained-active-only
  bbob-constrained-no-disguise)

## name | suite | observer | suite options
set(GOLDEN_OBSERVATIONS
  "bbob|bbob|bbob|dimensions: 2,5 function_indices: 1,8,15,24 instance_indices: 1,2"
  "bbob-old|bbob|bbob-old|dimensions: 2 function_indices: 1,8 instance_indices: 1"
  "bbob-boxed|bbob-boxed|bbob-boxed|dimensions: 2 function_indices: 1,15 instance_indices: 1"
  "bbob-noisy|bbob-noisy|bbob-noisy|dimensions: 2 function_indices: 1,30 instance_indices: 1"
  "bbob-largescale|bbob-largescale|bbob-largescale|dimensions: 20 function_indices: 1,10 instance_indices: 1"
  "bbob-mixint|bbob-mixint|bbob-mixint|dimensions: 5 function_indices: 1,15 instance_indices: 1"
  "bbob-biobj|bbob-biobj|bbob-biobj|dimensions: 2,5 function_indices: 1,20,55 instance_indices: 1,2"
  "bbob-biobj-ext|bbob-biobj-ext|bbob-biobj-ext|dimensions: 2 function_indices: 56,92 instance_indices: 1"
  "bbob-biobj-mixint|bbob-biobj-mixint|bbob-biobj-mixint|dimensions: 5 function_indices: 1,92 instance_indices: 1"
  "bbob-constrained|bbob-constrained|bbob-constrained|dimensions: 2,10 function_indices: 1,30,54 instance_indices: 1"
  "rw|bbob|rw|dimensions: 2 function_indices: 1 instance_indices: 1"
  "toy|toy|toy|dimensions: 2")

## Regenerate the references for this platform with
##   cmake --build build --target golden-update
## and review the changes before committing them.
add_custom_target(golden-update)

function(golden_scenario name)
  set(script_args
    -DGOLDEN=$<TARGET_FILE:golden>
    -DREFERENCE=${GOLDEN_REFERENCE_DIR}/${name}.txt
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/${name}
    ${ARGN})
  add_custom_target(golden-update-${name}
    COMMAND ${CMAKE_COMMAND} ${script_args} -DUPDATE=ON -P ${CMAKE_CURRENT_SOURCE_DIR}/golden.cmake
    DEPENDS golden)
  add_dependencies(golden-update golden-update-${name})
  if(EXISTS ${GOLDEN_REFERENCE_DIR})
    add_test(NAME golden-${name}
             COMMAND ${CMAKE_COMMAND} ${script_args} -P ${CMAKE_CURRENT_SOURCE_DIR}/golden.cmake)
  endif()
endfunction()

foreach(suite IN LISTS GOLDEN_SUITES)
  golden_scenario(values-${suite} -DMODE=values -DSUITE=${suite})
endforeach()

foreach(observation IN LISTS GOLDEN_OBSERVATIONS)
  string(REPLACE "|" ";" fields "${observation}")
  list(GET fields 0 name)
  list(GET fields 1 suite)
  list(GET fields 2 observer)
  list(GET fields 3 options)
  golden_scenario(observe-${name} -DMODE=observe -DSUITE=${suite} -DOBSERVER=${observer} "-DOPTIONS=${options}")
endforeach()

if(NOT EXISTS ${GOLDEN_REFERENCE_DIR})
  message(WARNING "No golden references for platform ${GOLDEN_PLATFORM}, golden tests are disabled. "
                  "Run the golden-update target to create them.")
endif()
