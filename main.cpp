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
#include "ConstraintGraph.h"
void print(const MemoryObject &m){
    if(m.getPtr()){
        llvm::outs() << *m.getPtr() << " (ID: " << m.getId() << ")";
    }
    else{
        llvm::outs() << "New MO (ID: " << m.getId() << ")";
    }
}


void printConstraints(const std::vector<Constraint> &constraints){
    for(auto cons : constraints){
        
        auto lhs = cons.getLhs();
        auto rhs = cons.getRhs();
        auto type = cons.getType();

        if(type == Constraint::ConstraintType::Load){
            // llvm::outs() << *rhs.getPtr() << " --- L ---> " << *lhs.getPtr() << "\n"; 
            print(lhs);
            llvm::outs() << " \t--- L --->\t ";
            print(rhs);
        }
        else if(type == Constraint::ConstraintType::Store){
            print(lhs);
            llvm::outs() << " \t--- S --->\t ";
            print(rhs);

        }
        else if(type == Constraint::ConstraintType::PointsTo){
            print(lhs);
            llvm::outs() << " \t--- P --->\t ";
            print(rhs);
        }
        else if(type == Constraint::ConstraintType::Copy){
            print(lhs);
            llvm::outs() << " \t--- C --->\t ";
            print(rhs);
        }

        llvm::outs() << "\n";
    }
}

MemoryObject getMemoryObjectFromPtr(std::vector<MemoryObject> &memoryObjects, const llvm::Value *ptr){
    auto iter = std::find_if(memoryObjects.begin(), memoryObjects.end(), [&ptr](const MemoryObject &mo) -> bool {return mo.getPtr() == ptr;});
    if(iter != memoryObjects.end()){
        return *iter;
    }
    else{
        auto mo = MemoryObject(ptr);
        memoryObjects.push_back(mo);
        return mo;
    }
}

ConstraintGraph createConstraintGraph(const std::vector<Constraint> &constraints){

    auto cg = ConstraintGraph(constraints);

    return cg;
}

int main(int argc, char** argv){

    llvm::LLVMContext Context;
    llvm::SMDiagnostic Err;

    // Read the LLVM IR file
    std::unique_ptr<llvm::Module> Mod = llvm::parseIRFile("/Users/jiaqi/Documents/PublicProject/GPU-FSPA/tests/misc/interSwap.ll", Err, Context);
    if (!Mod) {
        Err.print(argv[0], llvm::errs());
        return 1;
    }

    std::vector<MemoryObject> memoryObjects;
    std::vector<Constraint> constraints;

    // collect constraints
    for(const auto &Function : *Mod){
        for(const auto &BasicBlock : Function){
            for(const auto &Instruction : BasicBlock){
                // for x = alloca ..., introduce x points to m_x, where m_x is a new memoryobject allocated by this alloca (x, m_x, points-to)
                if(const auto &Alloca = llvm::dyn_cast<llvm::AllocaInst>(&Instruction)){
                    // create new memory object for lhs and rhs.
                    auto rhs = MemoryObject();
                    auto lhs = getMemoryObjectFromPtr(memoryObjects, Alloca);
                    memoryObjects.push_back(rhs);
                    // add constraints.
                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::PointsTo));

                }

                // for store x y, and x is a pointer, it means *y = x, add edge y --s--> x, (y, x, store)
                if(const auto &Store = llvm::dyn_cast<llvm::StoreInst>(&Instruction)){
                    const auto &PointerOperand = Store->getPointerOperand();
                    const auto &ValueOperand = Store->getValueOperand();

                    auto lhs = getMemoryObjectFromPtr(memoryObjects, PointerOperand);
                    auto rhs = getMemoryObjectFromPtr(memoryObjects, ValueOperand);
                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Store));
                }

                // for x = load y, it means x = *y, add edge y--l-->x, (y, x, load)
                if(const auto &Load = llvm::dyn_cast<llvm::LoadInst>(&Instruction)){
                    const auto &PointerOperand = Load->getPointerOperand();

                    auto lhs = getMemoryObjectFromPtr(memoryObjects, PointerOperand);
                    auto rhs = getMemoryObjectFromPtr(memoryObjects, Load);
                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Load));
                }
                
                // for x = bitcast y, it means x = y, add y--c-->x, (y, x, copy)
                if(const auto &BitCast = llvm::dyn_cast<llvm::BitCastInst>(&Instruction)){
                    const auto &PointerOperand = BitCast->getOperand(0);
                    if(!PointerOperand->getType()->isPointerTy()){
                        continue;
                    }

                    auto lhs = getMemoryObjectFromPtr(memoryObjects, PointerOperand);
                    auto rhs = getMemoryObjectFromPtr(memoryObjects, BitCast);
                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Copy));
                }
            
                if(const auto &Call = llvm::dyn_cast<llvm::CallInst>(&Instruction)){
                    // add copy edge from argument to parameter
                    // add copy edge from parameter to argument
                    size_t i = 0;
                    // llvm::outs() << "aaaa " << Call->getCalledFunction()->arg_size() << "\n";
                    while(Call->getCalledFunction() && i < Call->getCalledFunction()->arg_size()){

                        // todo: add case for indirect call
                        auto parameter = Call->getCalledFunction()->getArg(i);
                        auto argument = Call->getOperand(i);

                        // llvm::outs() << "aaaa" << *parameter << " " << *argument << "\n";

                        if(parameter->getType()->isPointerTy()){
                            auto lhs = getMemoryObjectFromPtr(memoryObjects, argument);
                            auto rhs = getMemoryObjectFromPtr(memoryObjects, parameter);
                            constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Copy));
                            constraints.push_back(Constraint(rhs, lhs, Constraint::ConstraintType::Copy));
                        }

                        ++i;
                    }
                    
                    
                }
            
            }
        }
    }


    printConstraints(constraints);
    auto cg = createConstraintGraph(constraints);

    llvm::outs() << "Node numbers: " << cg.getNodeNumbers() << ", p-edge number: " << cg.getPedgeNumbers() << ", c-edge number: " << cg.getCedgeNumbers()
        << ", s-edge number: " << cg.getSedgeNumbers() << ", l-edge number: " << cg.getLedgeNumbers() << "\n";

    // create worklist
    auto worklist = cg.getNodes();

    while(!worklist.empty()){
        std::set<MemoryObject> newWorklist;

        for(auto node : worklist){

            // print(node);
            // llvm::outs() << "\n";

            // apply store rule
            auto stores = std::set<MemoryObject>();
            if(cg.getSedges().count(node)){
                stores = cg.getSedges().at(node);
            }
            auto pts = std::set<MemoryObject>();
            if(cg.getPedges().count(node)){
                pts = cg.getPedges().at(node);
            }

            // for(auto store : stores){
            //     llvm::outs() << "store ";
            //     print(store);
            //     llvm::outs() << "\n";
            // }

            // for(auto pt : pts){
            //     llvm::outs() << "pt ";
            //     print(pt);
            //     llvm::outs() << "\n";
            // }




            for(auto store : stores){
                bool changed = false;
                for(auto pt : pts){
                    changed = cg.addCedge(store, pt);
                }
                if(changed){
                    newWorklist.insert(store);
                }

            }


            // apply load rule
            pts = std::set<MemoryObject>();
            if(cg.getPedges().count(node)){
                pts = cg.getPedges().at(node);
            }
            auto loads = std::set<MemoryObject>();
            if(cg.getLedges().count(node)){
                loads = cg.getLedges().at(node);
            }

            for(auto pt : pts){
                bool changed = false;
                for(auto load : loads){
                    changed = cg.addCedge(pt, load);
                }
                if(changed){
                    newWorklist.insert(pt);
                }
            }


            // apply copy rule
            auto copys = std::set<MemoryObject>();
            if(cg.getCedges().count(node)){
                copys = cg.getCedges().at(node);
            }
            pts = std::set<MemoryObject>();
            if(cg.getPedges().count(node)){
                pts = cg.getPedges().at(node);
            }

            for(auto copy : copys){
                bool changed = false;
                for(auto pt : pts){
                    changed = cg.addPedge(copy, pt);
                }
                if(changed){
                    newWorklist.insert(copy);
                }
            }

        }

        worklist = newWorklist;
    }

    // print detail of calculated cg
    for(auto p : cg.getPedges()){
        print(p.first);
        llvm::outs() << " ===>\n";

        for(auto pointee : p.second){
            llvm::outs() << "\t\t";
            print(pointee);
            llvm::outs() << "\n";
        }
    }

    // for(auto p : cg.getCedges()){
    //     print(p.first);
    //     llvm::outs() << " ===>\n";

    //     for(auto pointee : p.second){
    //         llvm::outs() << "\t\t";
    //         print(pointee);
    //         llvm::outs() << "\n";
    //     }
    // }

    // for(auto p : cg.getLedges()){
    //     print(p.first);
    //     llvm::outs() << " ===>\n";

    //     for(auto pointee : p.second){
    //         llvm::outs() << "\t\t";
    //         print(pointee);
    //         llvm::outs() << "\n";
    //     }
    // }


    return 0;
}