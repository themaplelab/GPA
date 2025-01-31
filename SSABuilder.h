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
    std::map<const llvm::BasicBlock*, std::set<const llvm::BasicBlock*>> DominanceFrontier;


    const llvm::BasicBlock *entry;
    std::set<const llvm::BasicBlock*> allBasicBlocks;

    public:

        // input of the constructor should be a cfg (or the entry block of a function since we are limiting ourselves within a function.)
        SSABuilder(const llvm::BasicBlock *entry) : entry(entry){

            findAllBasicBlocks();
            llvm::outs() << "Found all bbs:\n";
            for(auto p : bb2id){
                llvm::outs() << *(p.first->getFirstNonPHIOrDbg()) << " ==> " << p.second << "\n";
            }

            computeDominatorSet();

            computeImmediateDominator();

            buildDominatorTree();

            // DominatorTree = bb2TreeNode.at(entry); by far

            computeDominanceFrontier();
            llvm::outs() << "Dominance frontier:\n";
            for(auto p : DominanceFrontier){
                llvm::outs() << bb2id.at(p.first) << " ===>\n";
                for(auto d : p.second){
                    llvm::outs() << "\t" << bb2id.at(d) << "\n";;
                }
            }
        }

        void printDomTree();

    private:
        void buildDominatorTree();
        void findAllBasicBlocks();
        void computeDominatorSet();
        void computeImmediateDominator();
        void computeDominanceFrontier();
        void postOrderTraversal(std::unique_ptr<TreeNode>& root);

};




#endif