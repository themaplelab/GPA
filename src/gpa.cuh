#ifndef GPA_H
#define GPA_H

#include <string>



#ifdef __CUDACC__      // only visible when compiling as CUDA
__global__ void hello();
#endif

#ifdef __cplusplus     // safe for both CUDA and C++
extern "C" int gpamain(const std::string ptgFileName);
#endif






#endif