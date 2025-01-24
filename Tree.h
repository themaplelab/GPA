#ifndef GPA_TREE_H
#define GPA_TREE_H

#include <memory>
#include <vector>

#include "llvm/IR/BasicBlock.h"


class TreeNode{

    public:
        TreeNode(const llvm::BasicBlock *v) : val(v) {}
        void addChild(std::unique_ptr<TreeNode> &child){
            children.push_back(std::move(child));
        }

        std::vector<std::unique_ptr<TreeNode>>& getChildren() {return children;}
        const llvm::BasicBlock* getBasicBlock() {return val;}

    private:
        std::vector<std::unique_ptr<TreeNode>> children;
        const llvm::BasicBlock* val;

};






#endif