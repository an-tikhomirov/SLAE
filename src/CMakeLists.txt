set(HEADERS dense_matrix.hpp vector_operations.hpp csr_matrix.hpp FastGradientDescent.hpp SymGaussSeid.hpp SIMethods.hpp)
set(SOURCES dense_matrix.cpp csr_matrix.cpp vector_operations.cpp JacobiM.cpp SIM.cpp)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
