cmake_minimum_required(VERSION 3.19)
project(COCOIntegrationTest
        DESCRIPTION "COCO integration tests"
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
find_library(LIBRT rt)

if(MATH_LIBRARY)
    target_link_libraries(about_equal PUBLIC ${MATH_LIBRARY})
endif()

macro(integration_test name) 
  add_executable(${name} ${ARGN})  
  if (LIBM)
    target_link_libraries(${name} PUBLIC ${LIBM})
  endif()
  if (LIBRT)
    target_link_libraries(${name} PUBLIC ${LIBRT})
  endif()
  add_test(NAME ${name} 
           COMMAND $<TARGET_FILE:${name}>
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro()

integration_test(test_bbob-constrained test_bbob-constrained.c coco.cpp)
integration_test(test_bbob-largescale test_bbob-largescale.c coco.cpp)
integration_test(test_bbob-mixint test_bbob-mixint.cpp)
integration_test(test_biobj test_biobj.c coco.cpp)
integration_test(test_instance_extraction test_instance_extraction.c coco.cpp)