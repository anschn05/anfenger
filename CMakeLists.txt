cmake_minimum_required(VERSION 3.20)
project(ASC_bla)

set (CMAKE_CXX_STANDARD 17)

include_directories(src)

add_subdirectory (demos)


# Python bindings
if(SKBUILD)
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
  
  # Use pybind11 from pip install
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  find_package(pybind11 REQUIRED)
  
  pybind11_add_module(bla src/bind_bla.cpp)
  
  install(TARGETS bla DESTINATION ASCsoft)
endif()

# Fallback: allow manual build of Python module without scikit-build
if(NOT TARGET bla)
  find_package(Python COMPONENTS Interpreter Development REQUIRED)

    # Use pybind11 from pip install
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  
  find_package(pybind11 REQUIRED)
  pybind11_add_module(bla src/bind_bla.cpp)
endif()