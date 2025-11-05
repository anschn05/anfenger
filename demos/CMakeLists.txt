# demos/CMakeLists.txt
add_executable(demo_vector demo_vector.cpp)
add_executable(demo_matrix demo_matrix.cpp)

# Include-Pfade (damit "../src" nicht mehr nötig ist)
target_include_directories(demo_vector PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_include_directories(demo_matrix PUBLIC ${CMAKE_SOURCE_DIR}/src)

# Gemeinsame Ausgabe in demos/output/
set_target_properties(demo_vector demo_matrix PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/demos/output
)
