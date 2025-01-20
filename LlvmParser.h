#ifndef GPA_LLVM_PARSER_H
#define GPA_LLVM_PARSER_H

#include <memory>
#include <vector>
#include <set>
#include "llvm/IR/Module.h"
#include "llvm/IRReader/IRReader.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/Support/SourceMgr.h"

#include "llvm/Support/Signals.h"
#include "llvm/Support/PrettyStackTrace.h"
#include "Constraints.h"





class LLVMParser{

    std::unique_ptr<llvm::Module> module;
    std::set<size_t> topLevelVars;
    std::set<size_t> addressTakenVars;
    // All memory objects in this programs.
    std::vector<MemoryObject> memoryObjects;


    public:
        LLVMParser(std::string irName);
        std::unique_ptr<llvm::Module>& getLLVMModule() {return module;}
        std::vector<MemoryObject>& getMemoryObjects() {
            try{
                for(auto mo : memoryObjects){
                    llvm::outs() << mo.getPtr() << " " << *mo.getPtr() << " (ID: " << mo.getId() << ")\n";
                }
            }
            catch(...){
                llvm::sys::PrintStackTrace(llvm::errs());
            }
            
            llvm::outs() << "12312\n";
            return memoryObjects;
            }
        std::set<size_t> getTopLevelVariables() {return topLevelVars;}
        std::set<size_t> getAddressTakenVariables() {return addressTakenVars;}

    private:
        // implementation detail of public functions.

        void populateTopLevelAndAddressTakenVariables();
        size_t getMemoryObjectIndexFromPtr(const llvm::Value *, bool);

        

};



#endif