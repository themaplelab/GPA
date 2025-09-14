#ifndef GPA_SBV_GPU_H
#define GPA_SBV_GPU_H

#include <sys/types.h>
#include <iostream>
#include <climits>
#include <cassert>
#include <vector>
#include "Utils.cuh"



#define WORDS_PER_BLOCK 29
#define MAX_ALLOCATED_SBV_NUM 50000000	// 128 bytes SBV leads to about 8GB.
#define NULL_ADDR -1


using SBV_ADDR_TYPE = u_int64_t;

struct SBV;
__device__ u_int64_t gpu_malloc_sbv(SBV* pool, u_int32_t* pool_ptr);


struct SBV{

	// static SBV *pool;
	// static int *pool_ptr;

	int lock = 0;

	u_int32_t base;
	u_int32_t bits[WORDS_PER_BLOCK];
	u_int64_t next;

	__device__ __host__ void init(u_int32_t baseVal);

	__device__ __host__ bool isValidIndex(u_int32_t id);

	__device__ bool set(u_int32_t id, u_int64_t selfId, SBV *d_sbv_pool, u_int32_t *d_sbv_pool_ptr);

	__host__	bool setCPU(u_int32_t id, u_int64_t selfId, std::vector<SBV>& pool);

	__device__ bool test(u_int32_t id, u_int64_t selfId, SBV *d_sbv_pool);

	__device__ __host__	bool empty(SBV *d_sbv_pool) const;

	__device__ size_t toArray(uint32_t* out, SBV *d_sbv_pool, u_int32_t max_len);

	__host__ void print(const std::vector<SBV>& pool, u_int64_t headId) const;



};


struct SBVCPU{

	// static SBV *pool;
	// static int *pool_ptr;

	u_int32_t base;
	u_int32_t bits[WORDS_PER_BLOCK];
	SBVCPU *next;

	__device__ __host__ void init(u_int32_t baseVal);

	__device__ __host__ bool isValidIndex(u_int32_t id);

	// __device__ bool set(u_int32_t id, SBV *d_sbv_pool, u_int32_t *d_sbv_pool_ptr);

	__host__  bool setCPU(u_int32_t id);

	__device__ bool test(u_int32_t id);

	__device__ __host__ bool empty() const;

	__device__ int toArray(uint32_t* out, u_int32_t max_len);

	__host__ void print() const;



};




#endif 