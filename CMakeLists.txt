cmake_minimum_required(VERSION 3.20.0)

project(GPUFlowSensitivePointerAnalysis)


find_package(LLVM REQUIRED CONFIG)

message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")


include_directories(${LLVM_INCLUDE_DIRS})
separate_arguments(LLVM_DEFINITIONS_LIST NATIVE_COMMAND ${LLVM_DEFINITIONS})
add_definitions(${LLVM_DEFINITIONS_LIST})

add_executable(gpu-fspa main.cpp LlvmParser.cpp constraints.cpp AndersenPointerAnalysis.cpp SSABuilder.cpp)

llvm_map_components_to_libnames(llvm_libs support core irreader)

target_link_libraries(gpu-fspa ${llvm_libs})
target_compile_features(gpu-fspa PRIVATE cxx_std_17)
# target_compile_options(gpu-fspa PRIVATE -fsanitize=address)