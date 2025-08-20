#ifndef GPA_SBV_H
#define GPA_SBV_H

#include "Util/GeneralType.h"
#include "MSSA/MemSSA.h"
#include <iostream>
#include <set>

struct SBV{
	SVF::u32_t base;
	SVF::u32_t bits[29];
	SBV *next;

	SBV(SVF::u32_t b) : base(b){
		for(SVF::u32_t i = 0; i < 29; ++i){
			bits[i] = 0;
		}
	}
	SBV(){
		clear();
	}
	~SBV(){
		delete next;
	} 

	void clear(){
		base = 0;
		#pragma unroll
		for(SVF::u32_t i = 0; i < 29; ++i){
			bits[i] = 0;
		}
	}

	bool isValidIndex(SVF::u32_t id){
		return id > 0 && id < 928 * 0xFFFFFFFF;
	}

	bool set(SVF::u32_t id){

		bool changed = false;

		if(!isValidIndex(id)){
			return changed;
		}

		SVF::u32_t baseId = id / 928;
		SVF::u32_t bitsId = id % 928;
		SVF::u32_t wordId = bitsId / 32;
		SVF::u32_t inWordId = bitsId % 32;

		SBV *curr = this;
		SBV *prev = nullptr;

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
			SBV *newSbv = new SBV(baseId);
			newSbv->bits[28-wordId] |= (1u << inWordId);
			newSbv->next = curr;
			if(prev){
				prev->next = newSbv;
			}
			else{
				SBV *backup = new SBV(*this);

				this->base = baseId;
				for(SVF::u32_t i = 0; i < 29; ++i){
					this->bits[i] = 0;
				}
				this->set(id);
				this->next = backup;
			}
			changed = true;
		}
		return changed;
	}

	bool test(SVF::u32_t id){
		if(!isValidIndex(id)){
			return false;
		}

		SVF::u32_t baseId = id / 928;
		SVF::u32_t bitsId = id % 928;
		SVF::u32_t wordId = bitsId / 32;
		SVF::u32_t inWordId = bitsId % 32;

		const SBV *curr = this;
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

	bool unionWith(SBV *other){
		// merge two sorted linklist

		// todo: update unionwith to return bool

		SBV *s = other;
		SBV *curr = this;
		SBV *prev = nullptr;
		bool changed = false;

		while(s){
			while(curr && curr->base < s->base){
				prev = curr;
				curr = curr->next;
			}

			if(curr && curr->base == s->base){
				for(SVF::u32_t i = 0; i < 29; ++i){
					uint32_t old = curr->bits[i];
					curr->bits[i] |= s->bits[i];
					if(old != curr->bits[i]){
						changed = true;
					}
				}
				s = s->next;
			}
			else{
				if(prev){
					SBV *newSbv = new SBV(s->base);
					for(SVF::u32_t i = 0; i < 29; ++i){
						newSbv->bits[i] = s->bits[i];
					}
					newSbv->next = curr;

					prev->next = newSbv;
					prev = newSbv;
					s = s->next;
				}
				else{
					SBV *backup = new SBV(*this);
					this->base = s->base;
					for(SVF::u32_t i = 0; i < 29; ++i){
						this->bits[i] = s->bits[i];
					}
					this->next = backup;
					s = s->next;
				}
				changed = true;
			}
		}
		return changed;
	}

	void print() const{
		const SBV *curr = this;
		std::cout << "{";
		while(curr){
			SVF::u32_t id = curr->base * 928;
			for(SVF::u32_t i = 0; i < 29; ++i){
				if(curr->bits[28-i] == 0){
					continue;
				}
				for(SVF::u32_t j = 0; j < 32; ++j){
					if((curr->bits[28-i] >> j) & 1u){
						std::cout << id + i*32 + j << " ";
						// std::cout << i << " " << j << "\n";
					}
				}
			}

			curr = curr->next;
		}
		std::cout << "}\n";
	}

	void print(SVF::MemSSA *mssa) const{
		const SBV *curr = this;
		std::cout << "{";
		while(curr){
			SVF::u32_t id = curr->base * 928;
			for(SVF::u32_t i = 0; i < 29; ++i){
				if(curr->bits[28-i] == 0){
					continue;
				}
				for(SVF::u32_t j = 0; j < 32; ++j){
					if((curr->bits[28-i] >> j) & 1u){
						std::cout << mssa->getDetailedNameOfPtgNode(id + i*32 + j) << " ";
						// std::cout << i << " " << j << "\n";
					}
				}
			}

			curr = curr->next;
		}
		std::cout << "}\n";
	}

	void printRaw() const{
		const SBV *curr = this;
		while(curr){
			std::cout << curr->base << "\n";
			for(SVF::u32_t i = 0; i < 29; ++i){
				std::cout << curr->bits[i] << "\n";
			}
			curr = curr->next;
		}
	}

	std::set<SVF::u32_t> getNodeIds() const{
		std::set<SVF::u32_t> res;

		const SBV *curr = this;
		while(curr){
			SVF::u32_t id = curr->base * 928;
			for(SVF::u32_t i = 0; i < 29; ++i){
				if(curr->bits[28-i] == 0){
					continue;
				}
				for(SVF::u32_t j = 0; j < 32; ++j){
					if((curr->bits[28-i] >> j) & 1u){
						res.insert(id + i*32 + j);
					}
				}
			}
			curr = curr->next;
		}
		return res;

	}

	bool empty() const{
		for(SVF::u32_t i = 0; i < 29; ++i){
			if(bits[28-i] != 0){
				return false;
			}
		}
		return true;
	}

};




#endif