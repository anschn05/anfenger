add_executable(test_matrix_mul test_matrix_mul.cpp)
target_include_directories(test_matrix_mul PUBLIC ${CMAKE_SOURCE_DIR}/src)

# No extra libraries required; this is a small unit test for header-only code
