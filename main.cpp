#include "Constraints.h"
#include "ConstraintGraph.h"
#include "LlvmParser.h"
#include "PointerAnalysis.h"

#include <memory>

#include "llvm/IR/Module.h"
#include "llvm/Support/raw_ostream.h"


void print(const MemoryObject &m){
    if(m.getPtr()){
        llvm::outs() << *m.getPtr() << " (ID: " << m.getId() << ")";
    }
    else{
        llvm::outs() << "New MO (ID: " << m.getId() << ")";
    }
}

// todo: remove duplicate header
// todo: remove deadcode




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


    llvm::outs() << "Andersen pointer analysis\n";
    // collect constraints

    auto andersen = AndersenPointerAnalysis(parser, Mod);
    andersen.printConstraints();

    auto &cg = andersen.getConstraintGraph();

    llvm::outs() << "Node numbers: " << cg->getNodeNumbers() << ", p-edge number: " << cg->getPedgeNumbers() << ", c-edge number: " << cg->getCedgeNumbers()
        << ", s-edge number: " << cg->getSedgeNumbers() << ", l-edge number: " << cg->getLedgeNumbers() << "\n";

    andersen.analyze();
    
    // print detail of calculated cg
    andersen.printPEdge();
    andersen.printCEdge();
    andersen.printLEdge();
    andersen.printSEdge();



    return 0;
}



    

    


    

    



    


