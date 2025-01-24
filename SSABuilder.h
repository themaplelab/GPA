#ifndef GPA_SSA_BUILDER_H
#define GPA_SSA_BUILDER_H

#include "Tree.h"

#include <map>
#include <set>

#include "llvm/IR/BasicBlock.h"
#include "llvm/Support/raw_ostream.h"



class SSABuilder{

    std::map<const llvm::BasicBlock*, std::set<const llvm::BasicBlock*>> dominateSet;
    std::map<const llvm::BasicBlock*, const llvm::BasicBlock*> immediateDominatorMap;
    std::map<const llvm::BasicBlock*, std::set<const llvm::BasicBlock*>> immediateDominatedBy;
    std::map<const llvm::BasicBlock*, std::unique_ptr<TreeNode>> bb2TreeNode;
    std::map<const llvm::BasicBlock*, size_t> bb2id;

    const llvm::BasicBlock *entry;
    std::set<const llvm::BasicBlock*> allBasicBlocks;

    public:

        // input of the constructor should be a cfg (or the entry block of a function since we are limiting ourselves within a function.)
        SSABuilder(const llvm::BasicBlock *entry) : entry(entry){
            llvm::outs() << "1\n";
            findAllBasicBlocks();
            llvm::outs() << "2\n";

            computeDominatorSet();
            llvm::outs() << "3\n";

            computeImmediateDominator();
            llvm::outs() << "4\n";

            buildDominatorTree();
            llvm::outs() << "5\n";

            // DominatorTree = bb2TreeNode.at(entry); by far
        }

        void printDomTree();

    private:
        void buildDominatorTree();
        void findAllBasicBlocks();
        void computeDominatorSet();
        void computeImmediateDominator();

};




#endif