add_executable (demo_vector demo_vector.cpp)
target_sources (demo_vector PUBLIC ../src/vector.hpp)

# Add test_lapack executable which uses project headers and links to LAPACK/BLAS
add_executable(test_lapack test_lapack.cpp)
target_include_directories(test_lapack PUBLIC ${CMAKE_SOURCE_DIR}/src)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

if (DEFINED LAPACK_LIBRARIES)
	target_link_libraries(test_lapack PRIVATE ${LAPACK_LIBRARIES})
elseif (TARGET LAPACK::LAPACK)
	target_link_libraries(test_lapack PRIVATE LAPACK::LAPACK)
endif()

if (DEFINED BLAS_LIBRARIES)
	target_link_libraries(test_lapack PRIVATE ${BLAS_LIBRARIES})
elseif (TARGET BLAS::BLAS)
	target_link_libraries(test_lapack PRIVATE BLAS::BLAS)
endif()

# LU demo
add_executable(lu_demo lu_demo.cpp)
target_include_directories(lu_demo PUBLIC ${CMAKE_SOURCE_DIR}/src)
if (TARGET LAPACK::LAPACK)
	target_link_libraries(lu_demo PRIVATE LAPACK::LAPACK)
elseif (DEFINED LAPACK_LIBRARIES)
	target_link_libraries(lu_demo PRIVATE ${LAPACK_LIBRARIES})
endif()
if (TARGET BLAS::BLAS)
	target_link_libraries(lu_demo PRIVATE BLAS::BLAS)
elseif (DEFINED BLAS_LIBRARIES)
	target_link_libraries(lu_demo PRIVATE ${BLAS_LIBRARIES})
endif()
