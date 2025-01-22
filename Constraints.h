#ifndef GPU_FSPA_CONSTRAINTS_H
#define GPU_FSPA_CONSTRAINTS_H

#include <vector>

#include "llvm/IR/Value.h"
#include "llvm/Support/raw_ostream.h"


class MemoryObject{

    static size_t counter;

    size_t id;
    const llvm::Value *ptr;
    bool isAllocated;


    public:
        MemoryObject(const llvm::Value *ptr = nullptr, bool isAllocated = false) : ptr(ptr), isAllocated(isAllocated){
            // llvm::outs() << "constructor called\n";
            id = counter++;
        }

        const llvm::Value* getPtr() const{
            return ptr;
        }

        bool isAllocatedMemoryObject() const{
            return isAllocated;
        }

        size_t getId() const{
            return id;
        }

        bool operator==(const MemoryObject& other) const{
            return this->ptr == other.ptr && isAllocated == other.isAllocated;
        }

        static MemoryObject getSearchMemoryObject(const llvm::Value *ptr, bool isAllocated){
            return MemoryObject(0,ptr,isAllocated);
        }

        struct HashFunction{
            size_t operator()(const MemoryObject &mo) const
            {
            size_t xHash = std::hash<const void*>()(static_cast<const void*>(mo.ptr)) << 1;
            size_t yHash = std::hash<bool>()(mo.isAllocated);
            return xHash ^ yHash;
            }
        };

    private:
        MemoryObject(size_t id, const llvm::Value *ptr, bool isAllocated) : id(id), ptr(ptr), isAllocated(isAllocated) {}

};






class Constraint{

    public:
        enum ConstraintType{
            PointsTo, Copy, Store, Load, Define, Use, P_Define, P_Use, P_Copy, NonSSA, Already_Defined, Klique
        };

    private:
        size_t lhs;
        size_t rhs;
        ConstraintType type;

    public:
        Constraint(size_t lhs, size_t rhs, ConstraintType type) : lhs(lhs), rhs(rhs), type(type) {}
        size_t getLhs() const{
            return lhs;
        }

        size_t getRhs() const{
            return rhs;
        }

        ConstraintType getType() const{
            return type;
        }

};


#endif