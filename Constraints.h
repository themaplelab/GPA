#ifndef GPU_FSPA_CONSTRAINTS_H
#define GPU_FSPA_CONSTRAINTS_H

#include <vector>

#include "llvm/IR/Value.h"

class MemoryObject{

    static size_t counter;

    private:
        size_t id;
        const llvm::Value *ptr;


    public:
        MemoryObject(const llvm::Value *ptr = nullptr) : ptr(ptr){
            id = counter++;
        }

        const llvm::Value* getPtr() const{
            return ptr;
        }

        size_t getId() const{
            return id;
        }

        bool operator<(const MemoryObject& other) const{
            return this->id < other.id;
        }

};

size_t MemoryObject::counter = 0;




class Constraint{

    public:
        enum ConstraintType{
            PointsTo, Copy, Store, Load, Define, Use, P_Define, P_Use, P_Copy, NonSSA, Already_Defined, Klique
        };

    private:
        MemoryObject lhs;
        MemoryObject rhs;
        ConstraintType type;

    public:
        Constraint(MemoryObject lhs, MemoryObject rhs, ConstraintType type) : lhs(lhs), rhs(rhs), type(type) {}
        MemoryObject getLhs() const{
            return lhs;
        }

        MemoryObject getRhs() const{
            return rhs;
        }

        ConstraintType getType() const{
            return type;
        }

};


#endif