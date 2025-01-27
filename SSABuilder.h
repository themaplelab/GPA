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
            llvm::outs() << "Found all bbs:\n";
            for(auto p : bb2id){
                llvm::outs() << *(p.first->getFirstNonPHIOrDbg()) << " ==> " << p.second << "\n";
            }
            llvm::outs() << "2\n";

            computeDominatorSet();
            for(auto p : dominateSet){
                llvm::outs() << bb2id.at(p.first) << " dominated by \n";
                for(auto d : p.second){
                    llvm::outs() << "\t" << bb2id.at(d) << "\n";
                }
            }
            llvm::outs() << "3\n";

            computeImmediateDominator();
            for(auto p : immediateDominatorMap){
                llvm::outs() << bb2id.at(p.first) << " immediate dominated by \n";
                llvm::outs() << "\t" << bb2id.at(p.second) << "\n";
            }

            for(auto p : immediateDominatedBy){
                llvm::outs() << bb2id.at(p.first) << " immediate dominates \n";
                for(auto d : p.second){
                    llvm::outs() << "\t" << bb2id.at(d) << "\n";
                }
            }
            
            llvm::outs() << "4\n";

            buildDominatorTree();
            for(auto &p : bb2TreeNode){
                llvm::outs() << bb2id[p.first] << " " << bool(p.second) << "\n";
            }
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