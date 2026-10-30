cmake_minimum_required(VERSION 3.16)
project(ionff_terms LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(ionff_terms INTERFACE)
target_include_directories(ionff_terms INTERFACE include)

add_executable(ionff_eval src/ionff_eval.cpp)
target_link_libraries(ionff_eval PRIVATE ionff_terms)

enable_testing()
add_executable(test_terms tests/test_terms.cpp)
target_link_libraries(test_terms PRIVATE ionff_terms)
target_compile_options(test_terms PRIVATE -Wall -Wextra)
add_test(NAME terms_match_python COMMAND test_terms ${CMAKE_CURRENT_SOURCE_DIR}/tests/golden.csv)
