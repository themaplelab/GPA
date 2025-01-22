#ifndef GPA_LLVM_PARSER_H
#define GPA_LLVM_PARSER_H

#include "Constraints.h"

#include <memory>
#include <unordered_set>
#include <unordered_map>

#include "llvm/IR/Module.h"



class LLVMParser{
    llvm::LLVMContext context;
    std::unique_ptr<llvm::Module> module;
    
    std::unordered_set<size_t> topLevelVars;
    std::unordered_set<size_t> addressTakenVars;
    // All memory objects in this programs.
    std::unordered_set<MemoryObject, MemoryObject::HashFunction> memoryObjects;
    std::unordered_map<int, MemoryObject> id2MemoryObjects;


    public:
        LLVMParser(std::string irName);
        std::unique_ptr<llvm::Module>& getLLVMModule() {return module;}
        std::unordered_set<MemoryObject, MemoryObject::HashFunction>& getMemoryObjects() {return memoryObjects;}
        std::unordered_map<int, MemoryObject>& getMemoryObjectMap() {return id2MemoryObjects;}

        MemoryObject getMemoryObject(size_t id){
            return id2MemoryObjects[id];
        }
        std::unordered_set<size_t> getTopLevelVariables() {return topLevelVars;}
        std::unordered_set<size_t> getAddressTakenVariables() {return addressTakenVars;}
        size_t getMemoryObjectIndexFromPtr(const llvm::Value *, bool);

    private:
        // implementation detail of public functions.

        void populateTopLevelAndAddressTakenVariables();
        

        

};



#endif