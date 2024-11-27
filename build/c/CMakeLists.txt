cmake_minimum_required(VERSION 3.19)
project(coco
  DESCRIPTION "Example COCO experiment"
  LANGUAGES C CXX)
include(CTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_library(MATH_LIBRARY m)

## Add core coco experiment library to build
add_library(coco STATIC coco.cpp coco.h)
target_include_directories(coco PUBLIC .)
if(MATH_LIBRARY)
    target_link_libraries(coco PUBLIC ${MATH_LIBRARY})
endif()

## Add example experiment and link agains the coco library
## defined above.
add_executable(example_experiment example_experiment.c)
target_link_libraries(example_experiment PUBLIC coco)

## Add test harness 
add_executable(test_coco test_coco.c)
target_link_libraries(test_coco PUBLIC coco)

if(BUILD_TESTING)
  add_test(NAME bbob2009 
           COMMAND $<TARGET_FILE:test_coco> 
                   ${CMAKE_CURRENT_SOURCE_DIR}/bbob2009_testcases.txt)
endif()
