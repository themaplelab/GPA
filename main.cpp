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
#include "LlvmParser.h"
#include "AndersenPointerAnalysis.cpp"

void print(const MemoryObject &m){
    if(m.getPtr()){
        llvm::outs() << *m.getPtr() << " (ID: " << m.getId() << ")";
    }
    else{
        llvm::outs() << "New MO (ID: " << m.getId() << ")";
    }
}





int main(int argc, char** argv){

    auto parser = LLVMParser("/Users/jiaqi/Documents/PublicProject/GPU-FSPA/tests/misc/interSwap.ll");


    std::unique_ptr<llvm::Module> &Mod = parser.getLLVMModule();

    auto mos = parser.getMemoryObjects();

    llvm::outs() << "Top level vars:\n";
    for(auto tlv : parser.getTopLevelVariables()){
        llvm::outs() << tlv << " ";
        print(parser.getMemoryObject(tlv));
        llvm::outs() << "\n";
    }

    llvm::outs() << "Address taken vars:\n";
    for(auto &atv : parser.getAddressTakenVariables()){
        llvm::outs() << atv << " ";
        print(parser.getMemoryObject(atv));
        llvm::outs() << "\n";
    }

    // llvm::outs() << "Current memory objects in set:\n";
    // for(auto m : parser.getMemoryObjects()){
    //     llvm::outs() << *m.getPtr() << " (" << m.getId() << ", " << m.isAllocatedMemoryObject() << ")\n";
    // }

    // llvm::outs() << "Current memory objects in map:\n";
    // for(auto p : parser.getMemoryObjectMap()){
    //     llvm::outs() << p.first << " " << *(p.second).getPtr() << " (" << p.second.getId() << ", " << p.second.isAllocatedMemoryObject() << ")\n";
    // }

    // build inter-procedural call graph.


    // auto Andersen = AndersenPointerAnalysis(idug, mos);
    llvm::outs() << "Andersen pointer analysis\n";
    // collect constraints

    auto andersen = AndersenPointerAnalysis(parser, Mod);
    andersen.printConstraints();

    auto &cg = andersen.getConstraintGraph();

    llvm::outs() << "Node numbers: " << cg->getNodeNumbers() << ", p-edge number: " << cg->getPedgeNumbers() << ", c-edge number: " << cg->getCedgeNumbers()
        << ", s-edge number: " << cg->getSedgeNumbers() << ", l-edge number: " << cg->getLedgeNumbers() << "\n";
    



    return 0;
}



    

    


    

    // // create worklist
    // auto worklist = cg.getNodes();

    // while(!worklist.empty()){
    //     std::set<MemoryObject> newWorklist;

    //     for(auto node : worklist){

    //         // apply store rule
    //         auto stores = std::set<MemoryObject>();
    //         if(cg.getSedges().count(node)){
    //             stores = cg.getSedges().at(node);
    //         }
    //         auto pts = std::set<MemoryObject>();
    //         if(cg.getPedges().count(node)){
    //             pts = cg.getPedges().at(node);
    //         }

    //         for(auto store : stores){
    //             bool changed = false;
    //             for(auto pt : pts){
    //                 changed = cg.addCedge(store, pt);
    //             }
    //             if(changed){
    //                 newWorklist.insert(store);
    //             }

    //         }


    //         // apply load rule
    //         pts = std::set<MemoryObject>();
    //         if(cg.getPedges().count(node)){
    //             pts = cg.getPedges().at(node);
    //         }
    //         auto loads = std::set<MemoryObject>();
    //         if(cg.getLedges().count(node)){
    //             loads = cg.getLedges().at(node);
    //         }

    //         for(auto pt : pts){
    //             bool changed = false;
    //             for(auto load : loads){
    //                 changed = cg.addCedge(pt, load);
    //             }
    //             if(changed){
    //                 newWorklist.insert(pt);
    //             }
    //         }


    //         // apply copy rule
    //         auto copys = std::set<MemoryObject>();
    //         if(cg.getCedges().count(node)){
    //             copys = cg.getCedges().at(node);
    //         }
    //         pts = std::set<MemoryObject>();
    //         if(cg.getPedges().count(node)){
    //             pts = cg.getPedges().at(node);
    //         }

    //         for(auto copy : copys){
    //             bool changed = false;
    //             for(auto pt : pts){
    //                 changed = cg.addPedge(copy, pt);
    //             }
    //             if(changed){
    //                 newWorklist.insert(copy);
    //             }
    //         }

    //     }

    //     worklist = newWorklist;
    // }

    // print detail of calculated cg
    // for(auto p : cg.getPedges()){
    //     print(p.first);
    //     llvm::outs() << " ===>\n";

    //     for(auto pointee : p.second){
    //         llvm::outs() << "\t\t";
    //         print(pointee);
    //         llvm::outs() << "\n";
    //     }
    // }

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


