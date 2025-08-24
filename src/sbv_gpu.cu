#include "sbv_gpu.cuh"


__device__ __host__ void SBV::init(u_int32_t baseVal){
		base = baseVal;
		#pragma unroll
		for (int i = 0; i < WORDS_PER_BLOCK; ++i)
				bits[i] = 0;
		next = NULL_ADDR;
}

__device__ __host__ bool SBV::isValidIndex(u_int32_t id){
	return id > 0 && id < WORDS_PER_BLOCK * 32 * UINT32_MAX;
}

__device__ bool SBV::set(u_int32_t id, u_int64_t selfId, SBV *d_sbv_pool, u_int32_t *d_sbv_pool_ptr){
	
	if(!isValidIndex(id)){
		return false;
	}

	

	const u_int32_t baseId = id / 928;
	const u_int32_t bitsId = id % 928;
	const u_int32_t wordId = bitsId / 32;
	const u_int32_t inWordId = bitsId % 32;
	const u_int32_t mask = 1u << inWordId;

	u_int64_t curr = selfId;
	u_int64_t prev = NULL_ADDR;

	while(curr != NULL_ADDR && d_sbv_pool[curr].base < baseId){
		prev = curr;
		curr = d_sbv_pool[curr].next;
	}

	if(curr != NULL_ADDR && d_sbv_pool[curr].base == baseId){
		u_int32_t old = d_sbv_pool[curr].bits[28-wordId];
		atomicOr(&d_sbv_pool[curr].bits[28 - wordId], mask);
		// curr->bits[28-wordId] |= mask;
		return old != d_sbv_pool[curr].bits[28-wordId];
	}
	else{
		u_int64_t newIdx = gpu_malloc_sbv(d_sbv_pool, d_sbv_pool_ptr);
		SBV &newBlock = d_sbv_pool[newIdx];
		newBlock.base = baseId;
		#pragma unroll
		for(int i = 0; i < WORDS_PER_BLOCK; ++i){
			newBlock.bits[i] = 0;
		}
		// newBlock->bits[28 - wordId] |= mask;
		atomicOr(&newBlock.bits[28 - wordId], mask);
		newBlock.next = curr;

		if(prev != NULL_ADDR){
			d_sbv_pool[prev].next = newIdx;
		}
		else{
			u_int64_t backupIdx = gpu_malloc_sbv(d_sbv_pool, d_sbv_pool_ptr);
			SBV &backup = d_sbv_pool[backupIdx];
			SBV &head = d_sbv_pool[selfId];

			backup.base = head.base;
			#pragma unroll
			for(int i = 0; i < WORDS_PER_BLOCK; ++i){
				backup.bits[i] = head.bits[i];
			}
			backup.next = head.next;
			head.base = baseId;
			#pragma unroll
			for(int i = 0; i < WORDS_PER_BLOCK; ++i){
				head.bits[i] = 0;
			}
			// bits[28 - wordId] |= mask;
			atomicOr(&head.bits[28 - wordId], mask);
			head.next = backupIdx;
		}

		// unlockSBV(&this->lock);
		return true;
	}
}

__device__ bool SBV::test(u_int32_t id, u_int64_t selfId, SBV *d_sbv_pool){
	if(!isValidIndex(id)){
		return false;
	}

	u_int32_t baseId = id / 928;
	u_int32_t bitsId = id % 928;
	u_int32_t wordId = bitsId / 32;
	u_int32_t inWordId = bitsId % 32;

	u_int64_t curr = selfId;
	while(curr != NULL_ADDR){
		const SBV &node = d_sbv_pool[curr];
		if(node.base == baseId){
			return (node.bits[28-wordId] >> inWordId) & 1u;
		}
		if(node.base > baseId){
			break;
		}
		curr = node.next;
	}

	return false;
}

// __host__ __device__ bool emptyBlock() {
// 	#pragma unroll
// 	for (int i = 0; i < WORDS_PER_BLOCK; ++i) {
// 			// keep your original 28 - i order
// 			if (bits[28 - i] != 0u) return false;
// 	}
// 	return true;
// }

__device__ __host__ bool SBV::empty(SBV *d_sbv_pool) const{

	u_int64_t curr = static_cast<u_int64_t>(this - d_sbv_pool);

	while (curr != NULL_ADDR) {
		const SBV& node = d_sbv_pool[curr];
		for (u_int32_t i = 0; i < WORDS_PER_BLOCK; ++i) {
			if (node.bits[WORDS_PER_BLOCK - 1 - i] != 0) {
				return false;
			}
		}
		curr = node.next;
	}
	return true;

}

__device__ size_t SBV::toArray(uint32_t* out, SBV *d_sbv_pool, u_int32_t max_len = 4096){


	size_t count = 0;
	u_int64_t curr = static_cast<u_int64_t>(this - d_sbv_pool);
	

	while(curr != NULL_ADDR && count < max_len) {
		SBV &node = d_sbv_pool[curr];

		#pragma unroll
		for(int i = 0; i < WORDS_PER_BLOCK && count < max_len; ++i){
			uint32_t word = node.bits[WORDS_PER_BLOCK-1-i];
			if(word == 0){
				continue;
			}
			
			for(u_int32_t j = 0; j < 32 && count < max_len; ++j){
				if((word >> j) & 1u){
					out[count++] = node.base * 928 + i * 32 + j;
				}
			}
		}
		curr = node.next;
	}
	return count;  // number of valid elements in out[]
}

__host__	bool SBV::setCPU(u_int32_t id, u_int64_t selfId, std::vector<SBV>& pool){

	bool changed = false;

	if(!isValidIndex(id)){
		return changed;
	}

	u_int32_t baseId = id / 928;
	u_int32_t bitsId = id % 928;
	u_int32_t wordId = bitsId / 32;
	u_int32_t inWordId = bitsId % 32;


	u_int64_t curr = selfId;
	u_int64_t prev = NULL_ADDR;

	while(curr != NULL_ADDR && pool[curr].base < baseId){
		prev = curr;
		curr = pool[curr].next;
	}

	if(curr != NULL_ADDR && pool[curr].base == baseId){
		uint32_t old = pool[curr].bits[28-wordId];
		pool[curr].bits[28-wordId] |= (1u << inWordId);
		if(old != pool[curr].bits[28-wordId]){
			changed = true;
		}
	}
	else{
		// add new block
		const u_int64_t newIdx = static_cast<u_int64_t>(pool.size());
		SBV newSbv;
		newSbv.base = baseId;
		for(int i = 0; i < WORDS_PER_BLOCK; ++i){
			newSbv.bits[i] = 0u;
		} 
		newSbv.bits[28-wordId] |= (1u << inWordId);
		newSbv.next = curr;
		pool.emplace_back(newSbv);

		if(prev != NULL_ADDR){
			pool[prev].next = newIdx;
		}
		else{
			const u_int64_t backupIdx = static_cast<u_int64_t>(pool.size());
			pool.emplace_back(pool[selfId]);
			SBV &head = pool[selfId];
			head.base = baseId;
			for(int i = 0; i < WORDS_PER_BLOCK; ++i){
				head.bits[i] = 0;
			}
			head.bits[28-wordId] |= (1u << inWordId);
			head.next = backupIdx;
		}
		changed = true;
	}
	return changed;
}

__host__ void SBV::print(const std::vector<SBV>& pool, u_int64_t headId) const{
	
	u_int64_t curr = headId;
	std::cout << "{";
	while(curr != NULL_ADDR){
		const SBV &node = pool[curr];
		u_int32_t id = node.base * 928;
		for(u_int32_t i = 0; i < WORDS_PER_BLOCK; ++i){
			if(node.bits[28-i] == 0){
				continue;
			}
			for(u_int32_t j = 0; j < 32; ++j){
				if((node.bits[28-i] >> j) & 1u){
					std::cout << id + i*32 + j << " ";
				}
			}
		}
		curr = node.next;

	}
	std::cout << "}\n";
}





__device__ u_int64_t gpu_malloc_sbv(SBV* pool, u_int32_t* pool_ptr){
	assert(*pool_ptr < MAX_ALLOCATED_SBV_NUM && "Max number of node exceeded.");
	uint32_t idx = atomicAdd(pool_ptr, 1u);
  return static_cast<u_int64_t>(idx);
}



__device__ __host__ void SBVCPU::init(u_int32_t baseVal){
	base = baseVal;
	for(int i = 0; i < WORDS_PER_BLOCK; ++i)
		bits[i] = 0;
	next = nullptr;
}



// todo : check git log and copy the cpy version of sbv
__device__ __host__ bool SBVCPU::isValidIndex(u_int32_t id){
        return id > 0 && id < 928 * UINT32_MAX;
}

// __device__ bool SBVCPU::set(u_int32_t id, SBV *d_sbv_pool, u_int32_t *d_sbv_pool_ptr){

//         // lockSBV(&this->lock);

//         if(!isValidIndex(id)){
//                 return false;
//         }



//         const u_int32_t baseId = id / 928;
//         const u_int32_t bitsId = id % 928;
//         const u_int32_t wordId = bitsId / 32;
//         const u_int32_t inWordId = bitsId % 32;
//         const u_int32_t mask = 1u << inWordId;

//         SBVCPU *curr = this;
//         SBVCPU *prev = nullptr;

//         while(curr && curr->base < baseId){
//                 prev = curr;
//                 curr = curr->next;

//         }

//         if(curr && curr->base == baseId){
//                 u_int32_t old = curr->bits[28-wordId];
//                 atomicOr(&curr->bits[28 - wordId], mask);
//                 // curr->bits[28-wordId] |= mask;
//                 // unlockSBV(&this->lock);
//                 return old != curr->bits[28-wordId];
//         }
//         else{
//                 SBVCPU* newBlock = gpu_malloc_sbv(d_sbv_pool, d_sbv_pool_ptr);
//                 newBlock->base = baseId;
//                 for(int i = 0; i < WORDS_PER_BLOCK; ++i){
//                         newBlock->bits[i] = 0;
//                 }
//                 // newBlock->bits[28 - wordId] |= mask;
//                 atomicOr(&newBlock->bits[28 - wordId], mask);
//                 newBlock->next = curr;

//                 if(prev){
//                         prev->next = newBlock;
//                 }
//                 else{
//                         SBVCPU *backup = gpu_malloc_sbv(d_sbv_pool, d_sbv_pool_ptr);
//                         backup->base = base;
//                         for(int i = 0; i < WORDS_PER_BLOCK; ++i){
//                                 backup->bits[i] = bits[i];
//                         }
//                         backup->next = next;
//                         this->base = baseId;
//                         for(int i = 0; i < 29; ++i){
//                                 this->bits[i] = 0;
//                         }
//                         // bits[28 - wordId] |= mask;
//                         atomicOr(&bits[28 - wordId], mask);
//                         this->next = backup;
//                 }

//                 // unlockSBV(&this->lock);
//                 return true;
//         }
// }

__device__ bool SBVCPU::test(u_int32_t id){
        if(!isValidIndex(id)){
                return false;
        }

        u_int32_t baseId = id / 928;
        u_int32_t bitsId = id % 928;
        u_int32_t wordId = bitsId / 32;
        u_int32_t inWordId = bitsId % 32;

        const SBVCPU *curr = this;
        while(curr){
                if(curr->base == baseId){
                        return (curr->bits[28-wordId] >> inWordId) & 1u;
                }
                if(curr->base > baseId){
                        break;
                }
                curr = curr->next;
        }

        return false;
}

__device__ __host__ bool SBVCPU::empty() const{

        const SBVCPU* curr = this;
        while (curr) {
                for (u_int32_t i = 0; i < 29; ++i) {
                        if (curr->bits[28 - i] != 0) {
                                return false;
                        }
                }
                curr = curr->next;
        }
        return true;

}

__device__ int SBVCPU::toArray(uint32_t* out, u_int32_t max_len = 4096){
        int count = 0;
        SBVCPU* curr = this;

        while (curr) {
                for (int i = 0; i < WORDS_PER_BLOCK && count < max_len; ++i) {
                        if(curr->bits[28-i] == 0){
                                continue;
                        }
                        uint32_t word = curr->bits[28-i];
                        for(u_int32_t j = 0; j < 32 && count < max_len; ++j){
                                if((word >> j) & 1u){
                                        out[count++] = curr->base * 928 + i * 32 + j;
                                }
                        }
                }
                curr = curr->next;
        }
        return count;  // number of valid elements in out[]
}

__host__        bool SBVCPU::setCPU(u_int32_t id){

        bool changed = false;

        if(!isValidIndex(id)){
                return changed;
        }

        u_int32_t baseId = id / 928;
        u_int32_t bitsId = id % 928;
        u_int32_t wordId = bitsId / 32;
        u_int32_t inWordId = bitsId % 32;

        SBVCPU *curr = this;
        SBVCPU *prev = nullptr;

        while(curr && curr->base < baseId){
                prev = curr;
                curr = curr->next;
        }

        if(curr && curr->base == baseId){
                uint32_t old = curr->bits[28-wordId];
                curr->bits[28-wordId] |= (1u << inWordId);
                if(old != curr->bits[28-wordId]){
                        changed = true;
                }
        }
        else{
                // add new block
                SBVCPU *newSbv = new SBVCPU();
                newSbv->base = baseId;
                newSbv->bits[28-wordId] |= (1u << inWordId);
                newSbv->next = curr;
                if(prev){
                        prev->next = newSbv;
                }
                else{
                        SBVCPU *backup = new SBVCPU();
                        backup->base = base;
                        for(int i = 0; i < WORDS_PER_BLOCK; ++i){
                                backup->bits[i] = bits[i];
                        }
                        backup->next = next;

                        this->base = baseId;
                        for(int i = 0; i < 29; ++i){
                                this->bits[i] = 0;
                        }
                        this->setCPU(id);
                        this->next = backup;
                }
                changed = true;
        }
        return changed;
}

__host__ void SBVCPU::print() const{

        const SBVCPU *curr = this;
        std::cout << "{";
        while(curr){
                u_int32_t id = curr->base * 928;
                for(u_int32_t i = 0; i < 29; ++i){
                        if(curr->bits[28-i] == 0){
                                continue;
                        }
                        for(u_int32_t j = 0; j < 32; ++j){
                                if((curr->bits[28-i] >> j) & 1u){
                                        std::cout << id + i*32 + j << " ";
                                }
                        }
                }
                curr = curr->next;

        }
        std::cout << "}\n";
}