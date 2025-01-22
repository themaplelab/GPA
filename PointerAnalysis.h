#ifndef GPA_POINTER_ANALYSIS_H
#define GPA_POINTER_ANALYSIS_H

#include <vector>

#include "Constraints.h"
#include "LlvmParser.h"
#include <memory>
#include "llvm/IR/Module.h"
#include "ConstraintGraph.h"




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

    private:
        void collectConstraints(LLVMParser &parser, std::unique_ptr<llvm::Module> &Mod);
        void createConstraintGraph();
        
    
};





#endif