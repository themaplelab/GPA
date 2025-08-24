#ifndef GPA_H
#define GPA_H

#include <string>

struct PointsToGraph;

#ifdef __CUDACC__      // only visible when compiling as CUDA
__global__ void hello();
#endif

#ifdef __cplusplus     // safe for both CUDA and C++
extern "C" int gpamain(PointsToGraph &ptg, const std::string ptgFileName);
#endif






#endif