cmake_minimum_required(VERSION 3.19)
project(COCOUnitTest
        DESCRIPTION "COCO unit tests"
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

add_library(about_equal STATIC about_equal.c about_equal.h)
target_include_directories(about_equal PUBLIC .)
if(MATH_LIBRARY)
    target_link_libraries(about_equal PUBLIC ${MATH_LIBRARY})
endif()

macro(minunit_test name)  
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PUBLIC about_equal)
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

add_executable(test_tainted_problem test_tainted_problem.c coco.cpp)
if (LIBM)
  target_link_libraries(test_tainted_problem PUBLIC ${LIBM})
endif()
if (LIBRT)
  target_link_libraries(test_tainted_problem PUBLIC ${LIBRT})
endif()

## Tests that only use the public C interface
add_executable(test_coco_suite test_coco_suite.c coco.cpp)
if (LIBM)
  target_link_libraries(test_coco_suite PUBLIC ${LIBM})
endif()
add_test(NAME test_coco_suite
         COMMAND $<TARGET_FILE:test_coco_suite>
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME test_tainted-bbob COMMAND $<TARGET_FILE:test_tainted_problem> bbob)
set_property(TEST test_tainted-bbob PROPERTY WILL_FAIL YES)
add_test(NAME test_tainted-bbob-largescale COMMAND $<TARGET_FILE:test_tainted_problem> bbob-largescale)
set_property(TEST test_tainted-bbob-largescale PROPERTY WILL_FAIL YES)

minunit_test(test_biobj_utilities)
minunit_test(test_brentq)
minunit_test(test_coco_archive)
minunit_test(test_coco_observer)
minunit_test(test_coco_problem)
minunit_test(test_coco_string)
minunit_test(test_coco_utilities)
minunit_test(test_logger_bbob)
minunit_test(test_logger_biobj)
minunit_test(test_mo_utilities)