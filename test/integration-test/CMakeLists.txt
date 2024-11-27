cmake_minimum_required(VERSION 3.19)
project(COCOUnitTest
        DESCRIPTION "COCO CPP integration tests"
        LANGUAGES C CXX)

include(CTest)

find_library(LIBM m)
find_library(LIBRT rt)

if(MATH_LIBRARY)
    target_link_libraries(about_equal PUBLIC ${MATH_LIBRARY})
endif()

macro(integration_test name) 
  if(EXISTS "${PROJECT_SOURCE_DIR}/${name}.cpp")
    add_executable(${name} ${name}.cpp ${ARGN})
  else()
    add_executable(${name} ${name}.c ${ARGN})
  endif()
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

integration_test(test_bbob-constrained coco.cpp)
integration_test(test_bbob-largescale coco.cpp)
integration_test(test_bbob-mixint)
integration_test(test_biobj coco.cpp)
integration_test(test_instance_extraction coco.cpp)