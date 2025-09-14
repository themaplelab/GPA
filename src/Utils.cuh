#ifndef GPA_UTILS_H
#define GPA_UTILS_H

#include <iostream>

#define CUDA_CHECK(call)                                                         \
do {                                                                             \
    cudaError_t err = call;                                                     \
    if (err != cudaSuccess) {                                                   \
        std::cerr << "CUDA error in file " << __FILE__ << " at line : " <<			\
				__LINE__  << " " << cudaGetErrorString(err) << " \n";               		\
        exit(EXIT_FAILURE);                                                     \
    }                                                                            \
} while (0)


struct CudaTimer {
  cudaEvent_t start{}, stop{};
  explicit CudaTimer(unsigned flags = 0) {
    CUDA_CHECK(cudaEventCreateWithFlags(&start, flags));
    CUDA_CHECK(cudaEventCreateWithFlags(&stop,  flags));
    CUDA_CHECK(cudaEventRecord(start));
  }
  ~CudaTimer() {
    // best-effort cleanup; ignore errors in destructor
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  }
  // returns milliseconds; also synchronizes on stop
  float stop_ms() {
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float ms = 0.f;
    CUDA_CHECK(cudaEventElapsedTime(&ms, start, stop));
    return ms;
  }
};

template <class Fn>
float time_cuda_ms(Fn&& fn) {
  CudaTimer t;
  std::forward<Fn>(fn)();
  return t.stop_ms();
}



#endif