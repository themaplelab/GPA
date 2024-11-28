#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Module.h"
#include "llvm/IRReader/IRReader.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/Instructions.h"

#include <memory>
#include <vector>
#include <algorithm>

#include "Constraints.h"


int main(int argc, char** argv){

    llvm::LLVMContext Context;
    llvm::SMDiagnostic Err;

    // Read the LLVM IR file
    std::unique_ptr<llvm::Module> Mod = llvm::parseIRFile("/Users/jiaqi/Documents/PublicProject/GPU-FSPA/tests/misc/test.bc", Err, Context);
    if (!Mod) {
        Err.print(argv[0], llvm::errs());
        return 1;
    }

    std::vector<MemoryObject> memoryObjects;
    std::vector<Constraint> constraints;

    for(const auto &Function : *Mod){
        for(const auto &BasicBlock : Function){
            for(const auto &Instruction : BasicBlock){
                // y = alloca : y = &x

                if(const auto &Alloca = llvm::dyn_cast<llvm::AllocaInst>(&Instruction)){
                    // create new memory object for lhs and rhs.
                    auto rhs = MemoryObject();
                    auto lhs = MemoryObject(Alloca);

                    memoryObjects.push_back(lhs);
                    memoryObjects.push_back(rhs);
                    // add constraints.
                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::PointsTo));

                }

                // store x y : *y = x : y --s-> x

                if(const auto &Store = llvm::dyn_cast<llvm::StoreInst>(&Instruction)){
                    const auto &PointerOperand = Store->getPointerOperand();
                    const auto &ValueOperand = Store->getValueOperand();

                    auto iterLhs = std::find_if(memoryObjects.begin(), memoryObjects.end(), [&PointerOperand](const MemoryObject &mo) -> bool {return mo.getPtr() == PointerOperand;});
                    if(iterLhs != memoryObjects.end()){
                        // found
                        auto iterRhs = std::find_if(memoryObjects.begin(), memoryObjects.end(), [&ValueOperand](const MemoryObject &mo) -> bool {return mo.getPtr() == ValueOperand;});
                        if(iterRhs != memoryObjects.end()){
                            constraints.push_back(Constraint(*iterLhs, *iterRhs, Constraint::ConstraintType::Store));
                        }
                        else{
                            auto rhs = MemoryObject(ValueOperand);
                            constraints.push_back(Constraint(*iterLhs, rhs, Constraint::ConstraintType::Store));
                        }

                    }
                    else{

                        auto lhs = MemoryObject(PointerOperand);
                        auto iterRhs = std::find_if(memoryObjects.begin(), memoryObjects.end(), [&ValueOperand](const MemoryObject &mo) -> bool {return mo.getPtr() == ValueOperand;});
                        if(iterRhs != memoryObjects.end()){
                            constraints.push_back(Constraint(lhs, *iterRhs, Constraint::ConstraintType::Store));
                        }
                        else{
                            auto rhs = MemoryObject(Store);
                            constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Store));
                        }
                    }
                }
                // y = load x : y = *x
                // bitcast / gep : y = x
            }
        }
    }



    // Print the parsed module
    // Mod->print(llvm::outs(), nullptr);

    // std::cout << "Hello world.\n";

    return 0;
}