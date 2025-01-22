#ifndef GPA_POINTER_ANALYSIS_H
#define GPA_POINTER_ANALYSIS_H

#include "ConstraintGraph.h"
#include "Constraints.h"
#include "LlvmParser.h"

#include "llvm/IR/Module.h"

#include <vector>
#include <memory>





// Base class for all pointer analysis
class PointerAnalysis{

    public:


    private:


};

class AndersenPointerAnalysis : public PointerAnalysis{

    std::vector<Constraint> constraints;
    std::unique_ptr<ConstraintGraph> cg;

    public:
        AndersenPointerAnalysis(LLVMParser &parser, std::unique_ptr<llvm::Module> &Mod){
            collectConstraints(parser, Mod);
            createConstraintGraph();
        }
        void printConstraints();
        std::unique_ptr<ConstraintGraph>& getConstraintGraph(){
            return cg;
        }

        void analyze();
        void printPEdge();
        void printCEdge();
        void printLEdge();
        void printSEdge();


    private:
        void collectConstraints(LLVMParser &parser, std::unique_ptr<llvm::Module> &Mod);
        void createConstraintGraph();
        
    
};





#endif