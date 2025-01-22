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
        for(auto &arg : F.args()){
            // todo: only build parameter of pointer type.
            auto lhs = getMemoryObjectIndexFromPtr(&arg, false);
            topLevelVars.insert(lhs);
        }
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

    auto searchMo = MemoryObject::getSearchMemoryObject(val, isAllocated);
    auto iter = memoryObjects.find(searchMo);

    if(iter != memoryObjects.end()){
        return iter->getId();
    }
    else{
        auto mo = MemoryObject(val, isAllocated);
        memoryObjects.insert(mo);
        // using embrace is impoertant to avoid unintented creation of new memoryobject.
        id2MemoryObjects.try_emplace(mo.getId(), mo);

        return mo.getId();
    }
}
