#ifndef GPA_LLVM_PARSER_H
#define GPA_LLVM_PARSER_H

#include <memory>
#include <vector>
#include <set>
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
    std::vector<MemoryObject> memoryObjects;
    // todo: update modules to be a map. mapping each id to it's memoryobject (which contains the same id)
    std::vector<std::unique_ptr<llvm::Module>> modules;


    public:
        LLVMParser(std::string irName);
        std::unique_ptr<llvm::Module>& getLLVMModule() {return module;}
        std::vector<MemoryObject>& getMemoryObjects() {return memoryObjects;}
        std::set<size_t> getTopLevelVariables() {return topLevelVars;}
        std::set<size_t> getAddressTakenVariables() {return addressTakenVars;}
        // todo: update getMemoryObjectIndexFromPtr to return id instead of the position in a vector.
        size_t getMemoryObjectIndexFromPtr(const llvm::Value *, bool);

    private:
        // implementation detail of public functions.

        void populateTopLevelAndAddressTakenVariables();
        

        

};



#endif