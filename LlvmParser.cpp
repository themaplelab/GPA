#include "LlvmParser.h"


LLVMParser::LLVMParser(std::string irName){
    llvm::SMDiagnostic Err;
    module = llvm::parseIRFile(irName, Err, context);
    if (!module) {
        Err.print(irName.c_str(), llvm::errs());
        std::terminate();
    }

    populateTopLevelAndAddressTakenVariables();
}


void LLVMParser::populateTopLevelAndAddressTakenVariables(){
    for(auto &F : *module){
        for(auto &BB : F){
            for(auto &I : BB){
                if(const auto Alloca = llvm::dyn_cast<llvm::AllocaInst>(&I)){
                    // a = alloca type creates two memory objects.
                    auto lhs = getMemoryObjectIndexFromPtr(Alloca, false);
                    topLevelVars.insert(lhs);
                    auto rhs = getMemoryObjectIndexFromPtr(Alloca, true);
                    addressTakenVars.insert(rhs);

                }
                else if(I.getType()->isPointerTy()){
                    // we only interested in pointer related instructions.
                    auto lhs = getMemoryObjectIndexFromPtr(&I, false);
                    topLevelVars.insert(lhs);
                }
            }
        }
    }
}

/*
    return the unique id of the memoryobject corresponds to "val".
*/
size_t LLVMParser::getMemoryObjectIndexFromPtr(const llvm::Value *val, bool isAllocated){
    auto iter = std::find_if(memoryObjects.begin(), memoryObjects.end(), [&val, isAllocated](const MemoryObject &mo) -> bool {return mo.getPtr() == val && mo.isAllocatedMemoryObject() == isAllocated;});
    if(iter != memoryObjects.end()){
        return iter - memoryObjects.begin();
    }
    else{
        auto mo = MemoryObject(val, isAllocated);
        
        memoryObjects.push_back(mo);
        return memoryObjects.size()-1;
    }
}
