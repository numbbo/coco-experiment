cmake_minimum_required(VERSION 3.20)
 
project(COCO LANGUAGES C CXX)

## The COCO library is C++, the JNI glue code is C.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Disable floating-point contraction (FMA), see src/coco_fp_contract.h
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-ffp-contract=off COCO_HAVE_FP_CONTRACT_OFF)
if(COCO_HAVE_FP_CONTRACT_OFF)
  add_compile_options(-ffp-contract=off)
endif()

find_package(Java REQUIRED)
find_package(JNI REQUIRED)
include(UseJava)

add_jar(coco
  Benchmark.java 
  CocoJNI.java 
  ExampleExperiment.java
  Observer.java
  Problem.java
  Suite.java
  GENERATE_NATIVE_HEADERS coco-native)

add_library(CocoJNI SHARED CocoJNI.c coco.cpp)
target_link_libraries(CocoJNI PRIVATE coco-native)
target_include_directories(CocoJNI PRIVATE ${JNI_INCLUDE_DIRS})
