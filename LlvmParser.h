#ifndef GPA_LLVM_PARSER_H
#define GPA_LLVM_PARSER_H

#include <memory>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include "llvm/IR/Module.h"
#include "llvm/IRReader/IRReader.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/IR/Instructions.h"


#include "llvm/Support/Signals.h"
#include "llvm/Support/PrettyStackTrace.h"
#include "Constraints.h"





class LLVMParser{
    llvm::LLVMContext context;
    std::unique_ptr<llvm::Module> module;
    
    std::set<size_t> topLevelVars;
    std::set<size_t> addressTakenVars;
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
        std::set<size_t> getTopLevelVariables() {return topLevelVars;}
        std::set<size_t> getAddressTakenVariables() {return addressTakenVars;}
        size_t getMemoryObjectIndexFromPtr(const llvm::Value *, bool);

    private:
        // implementation detail of public functions.

        void populateTopLevelAndAddressTakenVariables();
        

        

};



#endif