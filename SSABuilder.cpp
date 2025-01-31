#include "SSABuilder.h"

#include <vector>
#include <queue>
#include <algorithm>
#include <iterator>

#include "llvm/IR/CFG.h"


void SSABuilder::printDomTree(){
    std::queue<const llvm::BasicBlock*> worklist;
    std::set<const llvm::BasicBlock*> visited;
    worklist.push(entry);
    llvm::outs() << "Dominator tree for function: " << entry->getParent()->getName().str() << "\n"; 

    while(!worklist.empty()){
        auto current = worklist.front();
        worklist.pop();
        auto &node = bb2TreeNode.at(current);
        


        if(visited.count(current)){
            continue;
        }
        visited.insert(current);

        for(auto child : node->getChildren()){
            llvm::outs() << bb2id.at(current) << " -- dom --> " << bb2id.at((child.get()->getBasicBlock())) << "\n";
            worklist.push(child.get()->getBasicBlock());
        }

    }


}



void SSABuilder::buildDominatorTree(){
    std::queue<const llvm::BasicBlock*> worklist;
    std::set<const llvm::BasicBlock*> visited;
    worklist.push(entry);

    while(!worklist.empty()){
        auto current = worklist.front();
        worklist.pop();

        if(visited.count(current)){
            continue;
        }
        visited.insert(current);

        for(auto child : immediateDominatedBy[current]){
            bb2TreeNode.at(current)->addChild(bb2TreeNode.at(child));

            worklist.push(child);
        }

    }
}

// todo : not tested
void SSABuilder::computeImmediateDominator(){
    for(auto BB : allBasicBlocks){
        // by this points, there should be no empty entry.
        auto dominators = dominateSet.at(BB);
        dominators.erase(BB);
        for(auto candidate : dominators){
            // if candidate only exists in dominateSet[candidate], it is the idom
            auto tmp = dominators;
            tmp.erase(candidate);
            bool exist = false;
            for(auto c : tmp){
                if(dominateSet[c].count(candidate)){
                    exist = true;
                    break;
                }
            }
            if(!exist){
                immediateDominatorMap[BB] = candidate;
                immediateDominatedBy[candidate].insert(BB);
                break;
            }
        }
    }
}

void SSABuilder::computeDominatorSet(){
    std::queue<const llvm::BasicBlock*> worklist;

    // initialization
    
    for(auto BB : allBasicBlocks){
        dominateSet.try_emplace(BB, std::set<const llvm::BasicBlock*>{BB});
        worklist.push(BB);

        auto node = std::make_unique<TreeNode>(BB);
        
        bb2TreeNode.try_emplace(BB, std::move(node));
        
    }

    while(!worklist.empty()){
        auto current = worklist.front();
        worklist.pop();
        // llvm::outs() << bb2id[current] << " has dom set \n";

        std::set<const llvm::BasicBlock*> dominatees;
        if(llvm::pred_begin(current) != llvm::pred_end(current)){
            dominatees = dominateSet.at(*(llvm::pred_begin(current)));
        }
        for(auto pred : llvm::predecessors(current)){
            std::set<const llvm::BasicBlock*> tmp;
            std::set_intersection(dominatees.begin(), dominatees.end(), dominateSet[pred].begin(), dominateSet[pred].end(), std::inserter(tmp, tmp.begin()));
            dominatees = tmp;
        }
        dominatees.insert(current);
        // for(auto d : dominatees){
        //     llvm::outs() << "\t" << bb2id[d] << "\n";
        // }

        auto oldSize = dominateSet[current].size();
        dominateSet[current] = dominatees;
        if(oldSize != dominateSet[current].size()){
            for(auto child : llvm::successors(current)){
                worklist.push(child);
            }
        }

    }


    
}


void SSABuilder::findAllBasicBlocks(){
    std::queue<const llvm::BasicBlock*> worklist;
    worklist.push(entry);
    size_t cnt = 0;

    while(!worklist.empty()){
        auto current = worklist.front();
        worklist.pop();
        if(allBasicBlocks.count(current)){
            continue;
        }
        allBasicBlocks.insert(current);
        bb2id.try_emplace(current, cnt++);

        for(auto child : llvm::successors(current)){
            worklist.push(child);
        }
    }
}

void SSABuilder::computeDominanceFrontier(){
    // Use dfs to achieve bottom-up traversal
    auto &root = bb2TreeNode.at(entry);
    postOrderTraversal(root);
}

void SSABuilder::postOrderTraversal(std::unique_ptr<TreeNode>& root){
    for(auto child : root->getChildren()){
        postOrderTraversal(child);
    }
    DominanceFrontier.try_emplace(root->getBasicBlock(), std::set<const llvm::BasicBlock*>{});
    for(auto node : llvm::successors(root->getBasicBlock())){
        if(immediateDominatorMap.at(node) != root->getBasicBlock()){
            DominanceFrontier.at(root->getBasicBlock()).insert(node);
        }
    }

    for(auto child : root->getChildren()){
        for(auto node : DominanceFrontier.at(child.get()->getBasicBlock())){
            if(immediateDominatorMap.at(node) != root->getBasicBlock()){
                DominanceFrontier.at(root->getBasicBlock()).insert(node);
            }
        }
    }
}