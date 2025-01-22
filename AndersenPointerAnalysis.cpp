#include "PointerAnalysis.h"


void AndersenPointerAnalysis::collectConstraints(LLVMParser &parser, std::unique_ptr<llvm::Module> &Mod){
    for(const auto &Function : *Mod){
        for(const auto &BasicBlock : Function){
            for(const auto &Instruction : BasicBlock){
                // for x = alloca ..., introduce x points to m_x, where m_x is a new memoryobject allocated by this alloca (x, m_x, points-to)
                if(const auto &Alloca = llvm::dyn_cast<llvm::AllocaInst>(&Instruction)){
                    // create new memory object for lhs and rhs.
                    auto rhs = parser.getMemoryObjectIndexFromPtr(Alloca, true);
                    auto lhs = parser.getMemoryObjectIndexFromPtr(Alloca, false);

                    llvm::outs() << Instruction << " generate constraint (" << lhs << ", " << rhs << ", p)\n"; 
                    // add constraints.
                    // todo: the id of a memoryobject is duplicated with the index returned by getMemoryObjectIndexFromPtr. remove either one.
                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::PointsTo));

                }

                // for store x y, and x is a pointer, it means *y = x, add edge y --s--> x, (y, x, store)
                if(const auto &Store = llvm::dyn_cast<llvm::StoreInst>(&Instruction)){
                    const auto &PointerOperand = Store->getPointerOperand();
                    const auto &ValueOperand = Store->getValueOperand();

                    // ignore storing non-pointer values for now.
                    if(!ValueOperand->getType()->isPointerTy()){
                        continue;
                    }

                    auto lhs = parser.getMemoryObjectIndexFromPtr(PointerOperand, false);
                    auto rhs = parser.getMemoryObjectIndexFromPtr(ValueOperand, false);
                    llvm::outs() << Instruction << " generate constraint (" << lhs << ", " << rhs << ", s)\n"; 

                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Store));
                }

                // for x = load y, it means x = *y, add edge y--l-->x, (y, x, load)
                if(const auto &Load = llvm::dyn_cast<llvm::LoadInst>(&Instruction)){
                    // ignore storing non-pointer values for now.
                    if(!Load->getType()->isPointerTy()){
                        continue;
                    }

                    const auto &PointerOperand = Load->getPointerOperand();

                    auto lhs = parser.getMemoryObjectIndexFromPtr(PointerOperand, false);
                    auto rhs = parser.getMemoryObjectIndexFromPtr(Load, false);
                    llvm::outs() << Instruction << " generate constraint (" << lhs << ", " << rhs << ", l)\n"; 

                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Load));
                }
                
                // for x = bitcast y, it means x = y, add y--c-->x, (y, x, copy)
                if(const auto &BitCast = llvm::dyn_cast<llvm::BitCastInst>(&Instruction)){
                    const auto &PointerOperand = BitCast->getOperand(0);
                    if(!PointerOperand->getType()->isPointerTy()){
                        continue;
                    }

                    auto lhs = parser.getMemoryObjectIndexFromPtr(PointerOperand, false);
                    auto rhs = parser.getMemoryObjectIndexFromPtr(BitCast, false);
                    llvm::outs() << Instruction << " generate constraint (" << lhs << ", " << rhs << ", c)\n"; 

                    
                    constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Copy));
                }
            
                if(const auto &Call = llvm::dyn_cast<llvm::CallInst>(&Instruction)){
                    // add copy edge from argument to parameter
                    // add copy edge from parameter to argument
                    size_t i = 0;
                    while(Call->getCalledFunction() && i < Call->getCalledFunction()->arg_size()){
                        auto parameter = Call->getCalledFunction()->getArg(i);
                        auto argument = Call->getOperand(i);

                        if(parameter->getType()->isPointerTy()){
                            auto lhs = parser.getMemoryObjectIndexFromPtr(argument, false);
                            auto rhs = parser.getMemoryObjectIndexFromPtr(parameter, false);
                            llvm::outs() << Instruction << " generate constraint (" << lhs << ", " << rhs << ", c)\n"; 

                            constraints.push_back(Constraint(lhs, rhs, Constraint::ConstraintType::Copy));
                            // constraints.push_back(Constraint(rhs, lhs, Constraint::ConstraintType::Copy));
                        }

                        ++i;
                    }
                }

                //todo : add return instructions handling
            
            }
        }
    }
}


void AndersenPointerAnalysis::printConstraints(){
    llvm::outs() << "Printing all constraints for Andersen pointer analysis\n";
    for(auto cons : constraints){
        
        auto lhs = cons.getLhs();
        auto rhs = cons.getRhs();
        auto type = cons.getType();

        if(type == Constraint::ConstraintType::Load){
            // llvm::outs() << *rhs.getPtr() << " --- L ---> " << *lhs.getPtr() << "\n"; 
            llvm::outs() << lhs << " \t--- L --->\t " << rhs;
        }
        else if(type == Constraint::ConstraintType::Store){
            llvm::outs() << lhs << " \t--- S --->\t " << rhs;

        }
        else if(type == Constraint::ConstraintType::PointsTo){
            llvm::outs() << lhs << " \t--- P --->\t " << rhs;
        }
        else if(type == Constraint::ConstraintType::Copy){
            llvm::outs() << lhs << " \t--- C --->\t " << rhs;
        }

        llvm::outs() << "\n";
    }
}


void AndersenPointerAnalysis::createConstraintGraph(){
    cg = std::make_unique<ConstraintGraph>(constraints);
}

void AndersenPointerAnalysis::analyze(){
    auto worklist = cg->getNodes();

        while(!worklist.empty()){
        std::set<size_t> newWorklist;

        for(auto node : worklist){

            // apply store rule
            auto stores = std::set<size_t>();
            if(cg->getSedges().count(node)){
                stores = cg->getSedges().at(node);
            }
            auto pts = std::set<size_t>();
            if(cg->getPedges().count(node)){
                pts = cg->getPedges().at(node);
            }

            for(auto store : stores){
                bool changed = false;
                for(auto pt : pts){
                    changed = cg->addCedge(store, pt);
                }
                if(changed){
                    newWorklist.insert(store);
                }

            }


            // apply load rule
            pts = std::set<size_t>();
            if(cg->getPedges().count(node)){
                pts = cg->getPedges().at(node);
            }
            auto loads = std::set<size_t>();
            if(cg->getLedges().count(node)){
                loads = cg->getLedges().at(node);
            }

            for(auto pt : pts){
                bool changed = false;
                for(auto load : loads){
                    changed = cg->addCedge(pt, load);
                }
                if(changed){
                    newWorklist.insert(pt);
                }
            }


            // apply copy rule
            auto copys = std::set<size_t>();
            if(cg->getCedges().count(node)){
                copys = cg->getCedges().at(node);
            }
            pts = std::set<size_t>();
            if(cg->getPedges().count(node)){
                pts = cg->getPedges().at(node);
            }

            for(auto copy : copys){
                bool changed = false;
                for(auto pt : pts){
                    changed = cg->addPedge(copy, pt);
                }
                if(changed){
                    newWorklist.insert(copy);
                }
            }

        }

        worklist = newWorklist;
    }
}

void AndersenPointerAnalysis::printPEdge(){
    for(auto p : cg->getPedges()){
        llvm::outs() << p.first << " ===P>\n";

        for(auto pointee : p.second){
            llvm::outs() << "\t\t" << pointee << "\n";;
        }
    }
}

void AndersenPointerAnalysis::printCEdge(){
    for(auto p : cg->getCedges()){
        llvm::outs() << p.first << " ===P>\n";

        for(auto pointee : p.second){
            llvm::outs() << "\t\t" << pointee << "\n";;
        }
    }
}

void AndersenPointerAnalysis::printLEdge(){
    for(auto p : cg->getLedges()){
        llvm::outs() << p.first << " ===P>\n";

        for(auto pointee : p.second){
            llvm::outs() << "\t\t" << pointee << "\n";;
        }
    }
}

void AndersenPointerAnalysis::printSEdge(){
    for(auto p : cg->getSedges()){
        llvm::outs() << p.first << " ===P>\n";

        for(auto pointee : p.second){
            llvm::outs() << "\t\t" << pointee << "\n";;
        }
    }
}