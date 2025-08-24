#include "Graphs/SVFG.h"
#include "SVF-LLVM/SVFIRBuilder.h"
#include "Util/Options.h"
#include "Util/SVFUtil.h"
#include "WPA/Andersen.h"
#include "MSSA/SVFGBuilder.h"
#include "sbv.h"
#include <map>
#include <string>
#include <chrono>
#include <fstream>
#include <filesystem>
#include "programl.h"
#include "SVFIR/SVFIR.h"
#include "gpa.cuh"



using namespace llvm;
using namespace std;
using namespace SVF;

/*!
 * An example to query alias results of two LLVM values
 */
SVF::AliasResult aliasQuery(PointerAnalysis *pta, const SVFVar *v1,
                            const SVFVar *v2) {
  return pta->alias(v1->getId(), v2->getId());
}

/*!
 * An example to print points-to set of an LLVM value
 */
std::string printPts(PointerAnalysis *pta, const SVFVar *svfval) {

  std::string str;
  raw_string_ostream rawstr(str);

  NodeID pNodeId = svfval->getId();
  const PointsTo &pts = pta->getPts(pNodeId);
  for (PointsTo::iterator ii = pts.begin(), ie = pts.end(); ii != ie; ii++) {
    rawstr << " " << *ii << " ";
    PAGNode *targetObj = pta->getPAG()->getGNode(*ii);
    rawstr << "(" << targetObj->toString() << ")\t ";
  }

  return rawstr.str();
}

/*!
 * An example to query/collect all successor nodes from a ICFGNode (iNode) along
 * control-flow graph (ICFG)
 */
void traverseOnICFG(ICFG *icfg, const Instruction *inst) {
  const ICFGNode *iNode = LLVMModuleSet::getLLVMModuleSet()->getICFGNode(inst);

  FIFOWorkList<const ICFGNode *> worklist;
  Set<const ICFGNode *> visited;
  worklist.push(iNode);

  /// Traverse along VFG
  while (!worklist.empty()) {
    const ICFGNode *iNode = worklist.pop();
    for (ICFGNode::const_iterator it = iNode->OutEdgeBegin(),
                                  eit = iNode->OutEdgeEnd();
         it != eit; ++it) {
      ICFGEdge *edge = *it;
      ICFGNode *succNode = edge->getDstNode();
      if (visited.find(succNode) == visited.end()) {
        visited.insert(succNode);
        worklist.push(succNode);
      }
    }
  }
}

void dummyVisit(const VFGNode *node) {}
/*!
 * An example to query/collect all the uses of a definition of a value along
 * value-flow graph (VFG)
 */
void traverseOnVFG(const SVFG *vfg, const SVFVar *svfval) {
  if (!vfg->hasDefSVFGNode(svfval))
    return;
  const VFGNode *vNode = vfg->getDefSVFGNode(svfval);
  FIFOWorkList<const VFGNode *> worklist;
  Set<const VFGNode *> visited;
  worklist.push(vNode);

  /// Traverse along VFG
  while (!worklist.empty()) {
    const VFGNode *vNode = worklist.pop();
    for (VFGNode::const_iterator it = vNode->OutEdgeBegin(),
                                 eit = vNode->OutEdgeEnd();
         it != eit; ++it) {
      VFGEdge *edge = *it;
      VFGNode *succNode = edge->getDstNode();
      if (visited.find(succNode) == visited.end()) {
        visited.insert(succNode);
        worklist.push(succNode);
      }
    }
  }

  /// Collect all LLVM Values
  for (Set<const VFGNode *>::const_iterator it = visited.begin(),
                                            eit = visited.end();
       it != eit; ++it) {
    const VFGNode *node = *it;
    dummyVisit(node);
    /// can only query VFGNode involving top-level pointers (starting with % or
    /// @ in LLVM IR) PAGNode* pNode = vfg->getLHSTopLevPtr(node); Value* val =
    /// pNode->getValue();
  }
}

void solve(SVF::MemSSA *mssa){
  // let's do this the dumb way
  // create a worklist
  // check each node in worklist check for possible rule application.

  outs() << "Performing fixed-point computation on set based points-to graph.\n";



  auto &ptg = mssa->getPtgMap();
  auto ptgNodeNumber = mssa->getPtgNodeNumber();

  // std::set<SVF::NodeID> worklist;

  // for(SVF::u32_t i = 0; i < ptgNodeNumber; ++i){
  //   worklist.insert(i);
  // }

  bool isChanged = true;

  auto start = std::chrono::high_resolution_clock::now();


  while(isChanged){

    isChanged = false;


    for(SVF::u32_t node = 0; node <= ptgNodeNumber; ++node){
          // update ptg to src => dst => type
      if(!ptg.count(node)){
        continue;
      }
      // outs() << "working on " << node << "\n";
      // rule 1
      // outs() << "rule1\n";
      if(ptg.at(node).count("c") && ptg.at(node).count("p")){
        for(auto cDst : ptg.at(node).at("c")){
          for(auto pDst : ptg.at(node).at("p")){
            if(mssa->addPtgEdge(cDst, pDst, "p")){
              isChanged = true;
            }
          }
        }
      }

      //rule 2
      // outs() << "rule2\n";
      if(ptg.at(node).count("u") && ptg.at(node).count("l")){
        for(auto uDst : ptg.at(node).at("u")){
          for(auto lDst : ptg.at(node).at("l")){
            if(ptg.count(uDst) && ptg.at(uDst).count("p_c")){
              for(auto pcDst : ptg.at(uDst).at("p_c")){
                if(pcDst == lDst){
                  if(mssa->addPtgEdge(uDst, pcDst, "c")){
                    isChanged = true;
                  }
                }
              }
            }
            
          }
        }
      }

      // rule 3
      // outs() << "rule3\n";
      if(ptg.at(node).count("s") && ptg.at(node).count("d")){
        // outs() << mssa->getDetailedNameOfPtgNode(node) << " has s and d edges\n";
        for(auto sDst : ptg.at(node).at("s")){
          for(auto dDst : ptg.at(node).at("d")){
            if(ptg.count(sDst) && ptg.at(sDst).count("p_c")){
              for(auto pcDst : ptg.at(sDst).at("p_c")){
                if(pcDst == dDst){
                  if(mssa->addPtgEdge(sDst, pcDst, "c")){
                    isChanged = true;
                  }
                }
              }
            }
          }
        }
      }

      // rule 4
      // outs() << "rule4\n";

      if(ptg.at(node).count("p_d") && ptg.at(node).count("p")){
        for(auto pdDst : ptg.at(node).at("p_d")){
          for(auto pDst : ptg.at(node).at("p")){
            if(mssa->isaVersionofPtr(pDst, pdDst)){
              if(mssa->addPtgEdge(node, pdDst, "d")){
                isChanged = true;
              }
            }
          }
        }
      }

      // rule 5
      // outs() << "rule5\n";
      if(ptg.at(node).count("p_u") && ptg.at(node).count("p")){
        for(auto puDst : ptg.at(node).at("p_u")){
          for(auto pDst : ptg.at(node).at("p")){
            if(mssa->isaVersionofPtr(pDst, puDst)){
              if(mssa->addPtgEdge(node, puDst, "u")){
                isChanged = true;
              }
            }
          }
        }
      }

      // rule 6
      // outs() << "rule6\n";

      if(ptg.at(node).count("d")){
        for(auto dDst : ptg.at(node).at("d")){
          if(ptg.count(dDst) && ptg.at(dDst).count("k")){
            for(auto kDst : ptg.at(dDst).at("k")){
              if(mssa->addPtgEdge(kDst, dDst, "a")){
                // should put all versions with n edges that end with dDst into worklist 
                isChanged = true;
              }
            }
          }
        }
      }

      // rule 7
      // outs() << "rule7\n";

      if(ptg.at(node).count("n")){
        for(auto nDst : ptg.at(node).at("n")){
          if(ptg.count(nDst) && ptg.at(nDst).count("k")){
            for(auto kDst : ptg.at(nDst).at("k")){
              if(ptg.count(kDst) && ptg.at(kDst).count("a")){
                for(auto aDst : ptg.at(kDst).at("a")){
                  if(mssa->getMrIdOfPtgNode(aDst) != mssa->getMrIdOfPtgNode(nDst)){
                    if(mssa->addPtgEdge(node, nDst, "c")){
                      isChanged = true;
                    }
                  }
                }
              }
            }
          }
        }
      }


    }


  }
  

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Runtime: " << duration.count() << " ms\n";
}

void sbvUnitTest(){
  // create sbv as points-to set
  // 1 -> {10, 20}
  std::map<SVF::NodeID, SBV*> pts;

  pts[1] = new SBV();
  // print empty
  pts[1]->print();


  pts[1]->set(10);
  pts[1]->set(20);
  // print 10 20
  pts[1]->print();
  // pts[1]->printRaw();

  // expected: 1 0 1
  outs() << pts[1]->test(10) << " " << pts[1]->test(11) << " " << pts[1]->test(20) << "\n"; 


  pts[0] = new SBV();
  pts[0]->unionWith(pts[1]);
  // expected 10 20
  pts[0]->print();

  

  // create another points-to set
  pts[2] = new SBV();
  pts[2]->set(10);
  // pts[2]->print();
  pts[2]->set(2000);
  // pts[2]->print();
  // pts[2]->printRaw();
  pts[2]->set(10000);
  // print 10 2000 10000
  pts[2]->print();


  // union two points-to sets
  pts[1]->unionWith(pts[2]);
  // print pts[1], expected 10 20 2000 10000
  pts[1]->print();
  // print pts[2], expected 10 2000 10000
  pts[2]->print();

  // test 3
  pts[3] = new SBV();
  pts[3]->set(30);
  pts[3]->set(2001);
  pts[3]->set(10001);
  pts[3]->set(9999);
  pts[2]->unionWith(pts[3]);

  // expected 10 20 30 2000 2001 9999 10000 10001
  pts[2]->print();
  // expected 30 2001 9999 10001
  pts[3]->print();


  return;
}

std::map<std::string, std::map<SVF::NodeID, SBV*>> translatePtgIntoSbv(SVF::MemSSA *mssa){
    std::map<std::string, std::map<SVF::NodeID, SBV*>> ptgSbv;

    auto &ptg = mssa->getPtgMap();
    auto ptgNodeNumber = mssa->getPtgNodeNumber();

    for(SVF::u32_t node = 0; node <= ptgNodeNumber; ++node){

      if(!ptg.count(node)){
        continue;
      }

      auto from = node;
      
      std::vector<std::string> edgeType{"p", "c", "s", "l", "p_d", "p_u", "p_c", "d", "u", "n", "k", "a"};
      for(auto et : edgeType){
        if(ptg.at(node).count(et)){
          for(auto to : ptg.at(node).at(et)){
            // std::cout << node << " " << et << " " << to << "\n";
            if(!ptgSbv[et].count(from)){
              ptgSbv[et][from] = new SBV();
            }
            
            ptgSbv[et][from]->set(to);
            // std::cout << "done\n";
          }
        }
      }
    }

    return ptgSbv;
}

void printSbvPts(std::map<std::string, std::map<SVF::NodeID, SBV*>> &ptgSbv, SVF::MemSSA *mssa){
  outs() << "printing sbv based pts\n";
  if(ptgSbv.count("p")){
    for(auto pair : ptgSbv.at("p")){
      if(pair.second->empty()){
        continue;
      }
      std::cout << mssa->getDetailedNameOfPtgNode(pair.first) << " ==> ";
      pair.second->print(mssa);
    }
  }
}

void solveSbv(std::map<std::string, std::map<SVF::NodeID, SBV*>> &ptgSbv, SVF::MemSSA *mssa){

  outs() << "Performing fixed-point computation for SBV based points-to graph.\n";
  bool changed = true;

  auto start = std::chrono::high_resolution_clock::now();

  while(changed){

    changed = false;


    // rule 1
    // outs() << "rule 1\n";
    if(ptgSbv.count("c") && ptgSbv.count("p")){
      for(auto pair1 : ptgSbv.at("c")){
        auto from1 = pair1.first;
        auto sbv1 = pair1.second;
        if(!ptgSbv.at("p").count(from1)){
          ptgSbv.at("p")[from1] = new SBV();
        }
        auto sbv2 = ptgSbv.at("p")[from1];
        for(auto nid : sbv1->getNodeIds()){
          if(!ptgSbv.count("p") || !ptgSbv.at("p").count(nid)){
            ptgSbv["p"][nid] = new SBV();
          }
          for(auto to : sbv2->getNodeIds()){
            if(ptgSbv["p"][nid]->set(to)){
              // std::cout << "Add p-edge (" << nid << "," << to << ")\n";
              changed = true;
            }
          }
        }
      }
    }

    // rule2
    // outs() << "rule 2\n";
    if(ptgSbv.count("u") && ptgSbv.count("l")){
      for(auto pair1 : ptgSbv.at("u")){
        auto from1 = pair1.first;
        for(auto pair2 : ptgSbv.at("l")){
          auto from2 = pair2.first;
          if(from1 == from2){
            auto sbv = pair1.second;
            for(auto nid : sbv->getNodeIds()){
              if(!ptgSbv.count("c") || !ptgSbv["c"].count(nid)){
                ptgSbv["c"][nid] = new SBV();
              }
              for(auto pcDst : ptgSbv["p_c"][nid]->getNodeIds()){
                if(ptgSbv.at("l")[from1]->test(pcDst)){
                  if(ptgSbv["c"][nid]->set(pcDst)){
                    // std::cout << "Add c-edge (" << nid << "," << pcDst << ")\n";
                    changed = true;
                  }
                }
              }
            }
          }
        }
      }
    }


    // rule 3
    // outs() << "rule 3\n";
    if(ptgSbv.count("s") && ptgSbv.count("d")){
      for(auto pair1 : ptgSbv.at("s")){
        auto from1 = pair1.first;
        auto sbv1 = pair1.second;
        if(ptgSbv.at("d").count(from1)){
          auto sbv2 = ptgSbv.at("d")[from1];
          for(auto sDst : sbv1->getNodeIds()){
            for(auto dDst : sbv2->getNodeIds()){
              if(ptgSbv.count("p_c") && ptgSbv.at("p_c").count(sDst)){
                if(ptgSbv.at("p_c").at(sDst)->test(dDst)){
                  if(!ptgSbv.count("c") || !ptgSbv.at("c").count(sDst)){
                    ptgSbv["c"][sDst] = new SBV();
                  }
                  if(ptgSbv.at("c").at(sDst)->set(dDst)){
                    // std::cout << "Add c-edge (" << sDst << "," << dDst << ")\n";
                    changed = true;
                  }
                }
              }
            }
          }
        }
        
      }
    }

    // rule 4
    // outs() << "rule 4\n";
    if(ptgSbv.count("p_d") && ptgSbv.count("p")){
      for(auto pair1 : ptgSbv.at("p_d")){
        auto from1 = pair1.first;
        auto sbv1 = pair1.second;
        if(ptgSbv.at("p").count(from1)){
          auto sbv2 = ptgSbv.at("p").at(from1);
          for(auto pdDst : sbv1->getNodeIds()){
            for(auto pDst : sbv2->getNodeIds()){
              if(mssa->isaVersionofPtr(pDst, pdDst)){
                if(!ptgSbv.count("d") || !ptgSbv["d"].count(from1)){
                  ptgSbv["d"][from1] = new SBV();
                }
                if(ptgSbv["d"][from1]->set(pdDst)){
                  // std::cout << "Add d-edge (" << from1 << "," << pdDst << ")\n";
                  changed = true;
                }
              }
            }
          }
        }
      }
    }

    // rule 5
    // outs() << "rule 5\n";
    if(ptgSbv.count("p_u") && ptgSbv.count("p")){
      for(auto pair1 : ptgSbv.at("p_u")){
        auto from1 = pair1.first;
        for(auto pair2 : ptgSbv.at("p")){
          auto from2 = pair2.first;
          if(from1 == from2){
            for(auto puDst : pair1.second->getNodeIds()){
              for(auto pDst : pair2.second->getNodeIds()){
                if(mssa->isaVersionofPtr(pDst, puDst)){
                  if(!ptgSbv.count("u") || !ptgSbv["u"].count(from1)){
                    ptgSbv["u"][from1] = new SBV();
                  }
                  if(ptgSbv["u"][from1]->set(puDst)){
                    // std::cout << "Add u-edge (" << from1 << "," << puDst << ")\n";
                    changed = true;
                  }
                }
              }
            }
          }
        }
      }
    }

    // rule 6
    // outs() << "rule 6\n";
    if(ptgSbv.count("d")){
      for(auto pair1 : ptgSbv.at("d")){
        for(auto dDst : pair1.second->getNodeIds()){
          if(ptgSbv.count("k") && ptgSbv["k"].count(dDst)){
            for(auto kDst : ptgSbv["k"][dDst]->getNodeIds()){
              if(!ptgSbv.count("a") || !ptgSbv["a"].count(kDst)){
                ptgSbv["a"][kDst] = new SBV();
              }
              if(ptgSbv["a"][kDst]->set(dDst)){
                // std::cout << "Add a-edge (" << kDst << "," << dDst << ")\n";
                changed = true;
              }
            }
          }
        }
      }
    }

    // rule 7
    // outs() << "rule 7\n";
    if(ptgSbv.count("n")){
      for(auto pair1 : ptgSbv.at("n")){
        auto from = pair1.first;
        for(auto nDst : pair1.second->getNodeIds()){
          if(ptgSbv.count("k") && ptgSbv["k"].count(nDst)){
            for(auto kDst : ptgSbv["k"][nDst]->getNodeIds()){
              if(ptgSbv.count("a") && ptgSbv["a"].count(kDst)){
                for(auto aDst : ptgSbv["a"][kDst]->getNodeIds()){
                  if(mssa->getMrIdOfPtgNode(aDst) != mssa->getMrIdOfPtgNode(nDst)){
                    if(!ptgSbv.count("c") || !ptgSbv["c"].count(from)){
                      ptgSbv["c"][from] = new SBV();
                    }
                    for(auto to : ptgSbv.at("n").at(from)->getNodeIds()){
                      if(ptgSbv["c"][from]->set(to)){
                        // std::cout << "Add c-edge (" << from << "," << to << ")\n";

                        changed = true;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Runtime: " << duration.count() << " ms\n";

}

void printPtgSbv(std::map<std::string, std::map<SVF::NodeID, SBV*>> &ptgSbv, SVF::MemSSA *mssa){
  outs() << "printing points-to set from SBV format\n";
  std::vector<std::string> edgeType{"p", "c", "s", "l", "p_d", "p_u", "p_c", "d", "u", "n", "k", "a"};
  for(auto et : edgeType){
    if(ptgSbv.count(et)){
      for(auto pair : ptgSbv.at(et)){
        if(pair.second->empty()){
          continue;
        }
        auto from = pair.first;
        std::cout << from << " == " << et << " ==> ";
        pair.second->print();
      }
    }
  }
}

void printRealPtg(SVF::MemSSA *mssa, bool detailed = true){
  // print points-to set
  auto realPtg = mssa->getPtgMap();
  outs() << "printing points-to set\n";
  for(SVF::NodeID i = 0; i < mssa->getPtgNodeNumber(); ++i){
    if(realPtg.count(i) && realPtg.at(i).count("p")){
      if(detailed){
        outs() << mssa->getDetailedNameOfPtgNode(i) << " ==> {";
      }
      else{
        outs() << i << " ==> {";
      }
      
      for(auto pointee : realPtg.at(i).at("p")){
        if(detailed){
        outs() << mssa->getDetailedNameOfPtgNode(pointee) << " ";
        }
        else{
          outs() << pointee << " ";
        }
        
      }
      outs() << "}\n";
    }
  }
}

void dumpRealPtg(SVF::MemSSA *mssa, std::string outputFileName){
  
  std::ofstream outFile(outputFileName+".out");
  if(!outFile){
    std::cerr << "Cannot open file " << outputFileName << ".out\n";
		std::terminate();
  }

  auto realPtg = mssa->getPtgMap();
  for(SVF::NodeID i = 0; i < mssa->getPtgNodeNumber(); ++i){
    if(realPtg.count(i) && realPtg.at(i).count("p")){
        outFile << i;
      for(auto pointee : realPtg.at(i).at("p")){
          outFile << "," << pointee;
      }
      outFile << "\n";
    }
  }

  outFile.close();
}

bool verifyPts(std::map<std::string, std::map<SVF::NodeID, SBV*>> &ptgSbv, SVF::MemSSA *mssa){
  // verify two computed pts are the same

  auto setBasedPtg = mssa->getTuplePtg();

  // traverse ptgSbv and remove entry in setBasedPtg. If not empty not entry not found, throw error.

  std::set<std::tuple<SVF::NodeID, SVF::NodeID, std::string>> extra;

  for(auto pair1 : ptgSbv){
    auto type = pair1.first;
    for(auto pair2 : pair1.second){
      auto from = pair2.first;

      for(auto to : pair2.second->getNodeIds()){
        std::tuple<SVF::NodeID, SVF::NodeID, std::string> t(from, to, type);
        auto iter = setBasedPtg.find(t);
        if(iter != setBasedPtg.end()){
          setBasedPtg.erase(iter);
        }
        else{
          extra.insert(t);
        }
      }
    }
  }

  bool error = false;
  if(!setBasedPtg.empty()){
    // outs() << "ptgSbv missing following entries:\n";
    // for(auto [from, to, type] : setBasedPtg){
    //   std::string fromName;
    //   std::stringstream ss(fromName); 
    //   ss << mssa->getDetailedNameOfPtgNode(from);
    //   ss << " == " << type << " ==> ";
    //   ss << mssa->getDetailedNameOfPtgNode(to);
    //   outs() << ss.str() << "\n";
    // }
    error = true;
  }

  if(!extra.empty()){
    // outs() << "ptgSbv contains following extra entries:\n";
    // for(auto [from, to, type] : extra){
    //   std::string fromName;
    //   std::stringstream ss(fromName); 
    //   ss << mssa->getDetailedNameOfPtgNode(from);
    //   ss << " == " << type << " ==> ";
    //   ss << mssa->getDetailedNameOfPtgNode(to);
    //   outs() << ss.str() << "\n";
    // }
    error = true;
  }

  if(error){
    outs() << "Fail to verify " << setBasedPtg.size() << " edges are missing, " <<  extra.size() << " extra edges\n";

    // printRealPtg(mssa);
    // printSbvPts(ptgSbv, mssa);
    return false;

  }
  else{
    outs() << "Verify succeed.\n";
  }

  return true;

}

void unitTest(){
  sbvUnitTest();
}

std::pair<std::vector<SVF::NodeID>, std::vector<SBV*>> getEdgeDataStructureForGpu(std::map<std::string, std::map<SVF::NodeID, SBV*>> &ptgSbv, std::string edgeType){
  std::vector<SVF::NodeID> nodeIds;
  std::vector<SBV*> sbvs;

  if(ptgSbv.count(edgeType)){
    for(auto pair : ptgSbv.at(edgeType)){
      nodeIds.push_back(pair.first);
      sbvs.push_back(pair.second);
    }

  }

  return std::make_pair(nodeIds, sbvs);

}

void dumpNodeId2mrIdAndVerMapToFile(std::string fileName, SVF::MemSSA *mssa){
  auto nodeId2mrIdAndVerMap = mssa->getNodeId2mrIdAndVerMap();
  std::ofstream nodeId2MrIdOutFile(fileName + "-nodeId2MrId.txt");
  if(!nodeId2MrIdOutFile){
    outs() << "ERROR! Cannot create file!\n";
    return;
  }

  for(auto pair : nodeId2mrIdAndVerMap){
    auto nodeId = pair.first;
    auto [from, to, type] = pair.second;
    nodeId2MrIdOutFile << nodeId << "," << from << "," << to << "," << type << "\n";
  }

  nodeId2MrIdOutFile.close();
}

void dumpNodeId2RelatedMrIdsMapToFile(std::string fileName, SVF::MemSSA *mssa){
  auto nodeId2RelatedMrIdsMap = mssa->getNodeId2MrIdsMap();
  std::ofstream nodeId2RelatedMrIdOutFile(fileName + "-nodeId2RelatedMrId.txt");
  if(!nodeId2RelatedMrIdOutFile){
    outs() << "ERROR! Cannot create file!\n";
    return;
  }

  for(auto pair : nodeId2RelatedMrIdsMap){
    auto nodeId = pair.first;
    nodeId2RelatedMrIdOutFile << nodeId;
    for(auto mrId : pair.second){
      nodeId2RelatedMrIdOutFile << "," << mrId;
    }
    nodeId2RelatedMrIdOutFile << "\n";
  }

  nodeId2RelatedMrIdOutFile.close();

}

void dumpPtgToFile(std::string fileName, SVF::MemSSA *mssa){
  auto ptg = mssa->getTuplePtg();

  std::ofstream outFile(fileName + ".txt");
  if(!outFile){
    outs() << "ERROR! Cannot create file!\n";
    return;
  }

  for(auto [from, to, type] : ptg){
    outFile << from << "," << to << "," << type << "\n";
  }

  outFile.close();
}

void dumpGroundTruthPtgToFile(std::string fileLoc, std::string fileName, SVF::MemSSA *mssa){
  auto ptg = mssa->getTuplePtg();

  std::filesystem::path p = fileLoc;
  auto parentDir = p.parent_path();
  parentDir /= "outputs";
  parentDir /= fileName;
  std::string fullPath = parentDir.string();

  // std::ofstream outFile(fullPath + ".out");
  std::ofstream outFile(fileName + ".out");

  if(!outFile){
    outs() << "ERROR! Cannot create file!\n";
    return;
  }

  for(auto [from, to, type] : ptg){
    outFile << from << "," << to << "," << type << "\n";
  }

  outFile.close();
}

void dumpNecessaryFilesForGPU(std::string fileLoc, std::string fileName, SVF::MemSSA *mssa){
  std::filesystem::path p = fileLoc;
  auto parentDir = p.parent_path();
  parentDir /= "outputs";
  parentDir /= fileName;
  std::string fullPath = parentDir.string();
  // std::cout << fullPath << "\n";

  // dumpPtgToFile(fullPath+"-initialPtg", mssa);
  // dumpNodeId2RelatedMrIdsMapToFile(fullPath, mssa);
  // dumpNodeId2mrIdAndVerMapToFile(fullPath, mssa);

  dumpPtgToFile(fileName+"-initialPtg", mssa);
  dumpNodeId2RelatedMrIdsMapToFile(fileName, mssa);
  dumpNodeId2mrIdAndVerMapToFile(fileName, mssa);
}

std::string getGraphOutputName(std::string filename){
  std::filesystem::path filePath = filename;
  return filePath.stem().string();
}

std::string getGraphOutputLocation(std::string filename){
  std::filesystem::path filePath = filename;
  return filePath.parent_path().string();
}

ProGraML generateExtendedProGraMlGraph(SVF::MemSSA *mssa, Andersen *ander, SVFIR *pag){

  ProGraML graph;
  std::vector<const SVFVar*> allVars;

  

  const CallGraph* svfirCallGraph = pag->getCallGraph();

  for(const auto& item: *svfirCallGraph){


    const FunObjVar* fun = item.second->getFunction();
    if(Options::MSSAFun()!="" && Options::MSSAFun()!=fun->getName()){
        continue;
    }
    if(fun->isDeclaration()){
      continue;
    }


    auto entry = graph.addGraphNodeFromIcfgNode(pag->getICFG()->getFunEntryICFGNode(fun));
    graph.addEdge(0, entry, "call");
    auto entryOutEdges = pag->getICFG()->getFunEntryICFGNode(fun)->getOutEdges();

    for(auto oe : entryOutEdges){
      auto toId = graph.addGraphNodeFromIcfgNode(oe->getDstNode());
      graph.addEdge(entry, toId, "control");
    }


    auto exit = graph.addGraphNodeFromIcfgNode(pag->getICFG()->getFunExitICFGNode(fun));
    graph.addEdge(exit, 0, "call");
    auto exitOutEdges = pag->getICFG()->getFunExitICFGNode(fun)->getOutEdges();
    for(auto oe : exitOutEdges){
      auto toId = graph.addGraphNodeFromIcfgNode(oe->getDstNode());
      graph.addEdge(exit, toId, "call");

      if(isa<RetICFGNode>(oe->getDstNode())){
        auto exitOutEdges = oe->getDstNode()->getOutEdges();
        for(auto subOe : exitOutEdges){
          auto subToId = graph.addGraphNodeFromIcfgNode(subOe->getDstNode());
          graph.addEdge(toId, subToId, "control");
        }
      }
    }
    

    if(mssa->hasFuncEntryChi(fun)){
      MemSSA::CHISet &entry_chis = mssa->getFuncEntryChiSet(fun);
      for(MemSSA::CHISet::iterator chi_it = entry_chis.begin(); chi_it != entry_chis.end(); chi_it++){
          auto chi = *chi_it;

          auto mrId = chi->getMR()->getMRID();
          auto resMrVer = chi->getResVer()->getSSAVersion();
          auto opMrVer = chi->getOpVer()->getSSAVersion();
          auto pts = chi->getMR()->getPointsTo();
          // for(NodeBS::iterator ii = pts.begin(), ie = pts.end(); ii != ie; ii++){
          //     nodeId2MrIdsMap[*ii].insert(mrId);
          // }

          auto mrDst = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,resMrVer);
          auto mrSrc = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,opMrVer);

          graph.addEdge(mrSrc, mrDst, "data");


          
      }
    } 


    for(FunObjVar::const_bb_iterator bit = fun->begin(), ebit = fun->end(); bit != ebit; ++bit){

      const SVFBasicBlock* bb = bit->second;

      std::vector<const SVFVar*> allVarInFun;

      MemSSA::PHISet& phiSet = mssa->getPHISet(bb);
      for(MemSSA::PHISet::iterator pi = phiSet.begin(), epi = phiSet.end(); pi !=epi; ++pi){

          // MemSSA::PHI phinode = *pi;

          auto mrId = (*pi)->getMR()->getMRID();
          auto resPiVer = (*pi)->getResVer()->getSSAVersion();
          auto resPiId = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,resPiVer);
          for(auto it = (*pi)->opVerBegin(), eit = (*pi)->opVerEnd(); it!=eit; ++it){
              auto opPiVer = it->second->getSSAVersion();
              auto opPiId = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,opPiVer);
              graph.addEdge(opPiId, resPiId, "data");
              // addPtgEdge(opPiId, resPiId, "c");
          }
      }


      for(const auto& inst: bb->getICFGNodeList()){
          bool isAppCall = SVFUtil::isNonInstricCallSite(inst) && !SVFUtil::isExtCall(inst);
          auto instId = graph.addGraphNodeFromIcfgNode(inst);
          if(isAppCall || SVFUtil::isHeapAllocExtCall(inst)){
            const CallICFGNode* cs = cast<CallICFGNode>(inst);
            

            for(auto stmt : inst->getSVFStmts()){
              if(auto call = dyn_cast<CallPE>(stmt)){
                auto src = call->getRHSVar();
                auto dst = call->getLHSVar();

                auto rhsId = graph.addGraphNodeFromSvfVar(src);
                auto lhsId = graph.addGraphNodeFromSvfVar(dst);
                // todo add edge label.
                graph.addEdge(rhsId, instId, "data");
                graph.addEdge(instId, lhsId, "data");
              }
            }

            auto fromId = graph.addGraphNodeFromIcfgNode(inst);
            auto outEdges = inst->getOutEdges();
            for(auto oe : outEdges){
              auto toId = graph.addGraphNodeFromIcfgNode(oe->getDstNode());
              graph.addEdge(fromId, toId, "call");
            }

            if(cs->isIndirectCall()){
              continue;
            }


            if(mssa->hasMU(cs)){
              for(MemSSA::MUSet::iterator mit = mssa->getMUSet(cs).begin(), emit = mssa->getMUSet(cs).end(); mit != emit; ++mit){

                auto mrId = (*mit)->getMR()->getMRID();
                auto mrVer = (*mit)->getMRVer()->getSSAVersion();

                auto mrSrcId = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,mrVer);
                auto mrDstId = graph.addGraphNodeFromMemoryRegion(cs->getCalledFunction()->getName(),mrId,0);

                graph.addEdge(mrSrcId, mrDstId, "data");
              }
            }

            if(mssa->hasCHI(cs)){
              for(MemSSA::CHISet::iterator cit = mssa->getCHISet(cs).begin(), ecit = mssa->getCHISet(cs).end(); cit != ecit; ++cit){
                auto mrId = (*cit)->getMR()->getMRID();
                auto mrResVer = (*cit)->getResVer()->getSSAVersion();

                auto mrDstId = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,mrResVer);

                auto mrOpVer = (*cit)->getOpVer()->getSSAVersion();
                auto mrSrcId = graph.addGraphNodeFromMemoryRegion(cs->getCalledFunction()->getName(),mrId,mrOpVer);

                graph.addEdge(mrSrcId, mrDstId, "def");
              }
            }
          }
          else{
            auto fromId = graph.addGraphNodeFromIcfgNode(inst);
            auto outEdges = inst->getOutEdges();
            for(auto oe : outEdges){
              auto toId = graph.addGraphNodeFromIcfgNode(oe->getDstNode());
              graph.addEdge(fromId, toId, "control");
            }

            // add data flow
            // for each inst, e.g. a = load b, get node for a, create data edge to a, get node for b, create edge from b.
            for(auto stmt : inst->getSVFStmts()){
              if(auto alloca = dyn_cast<AddrStmt>(stmt)){
                auto top = alloca->getLHSVar();
                auto topId = graph.addGraphNodeFromSvfVar(top);
                auto addr = alloca->getRHSVar();
                auto addrId = graph.addGraphNodeFromSvfVar(addr);
                graph.addEdge(instId, topId, "data");
                graph.addEdge(instId, addrId, "data");

                allVarInFun.push_back(top);
                allVars.push_back(top);

                auto m = mssa->getNodeId2MrIdsMap();
                if(m.count(alloca->getRHSVarID())){
                  for(auto mrId : m.at(alloca->getRHSVarID())){
                    auto mrDst = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,0);
                    graph.addEdge(addrId, mrDst, "data");
                  }
                }
                

              }
              else if(auto multiOp = dyn_cast<MultiOpndStmt>(stmt)){
                auto dst = multiOp->getRes();
                allVarInFun.push_back(dst);
                allVars.push_back(dst);

                auto dstId = graph.addGraphNodeFromSvfVar(dst);
                for(auto src : multiOp->getOpndVars()){
                  auto srcId = graph.addGraphNodeFromSvfVar(src);
                  graph.addEdge(srcId, instId, "data");
                }
                graph.addEdge(instId, dstId, "data");
              }
              else if(auto store = dyn_cast<StoreStmt>(stmt)){
                auto src = store->getRHSVar();
                auto dst = store->getLHSVar();

                auto rhsId = graph.addGraphNodeFromSvfVar(src);
                auto lhsId = graph.addGraphNodeFromSvfVar(dst);

                graph.addEdge(rhsId, instId, "data");
                graph.addEdge(lhsId, instId, "data");

                MemSSA::CHISet& chiSet = mssa->getCHISet(store);
                for(MemSSA::CHISet::iterator it = chiSet.begin(), eit = chiSet.end(); it!=eit; ++it){

                  auto mrId = (*it)->getMR()->getMRID();
                  auto mrOpVer = (*it)->getOpVer()->getSSAVersion();
                  auto mrResVer = (*it)->getResVer()->getSSAVersion();

                  auto mrOpId = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,mrOpVer);
                  auto mrResId = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,mrResVer);

                  graph.addEdge(mrOpId, instId, "def");
                  graph.addEdge(instId, mrResId, "def");
                }
              }
              else if(auto load = dyn_cast<LoadStmt>(stmt)){
                auto src = load->getRHSVar();
                auto dst = load->getLHSVar();
                allVarInFun.push_back(dst);
                allVars.push_back(dst);


                auto rhsId = graph.addGraphNodeFromSvfVar(src);
                auto lhsId = graph.addGraphNodeFromSvfVar(dst);
                // todo add edge label.
                graph.addEdge(rhsId, instId, "data");
                graph.addEdge(instId, lhsId, "data");

                MemSSA::MUSet &muSet = mssa->getMUSet(load);

                for(MemSSA::MUSet::iterator it = muSet.begin(), eit = muSet.end(); it!=eit; ++it){

                  auto mrId = (*it)->getMR()->getMRID();
                  auto mrVer = (*it)->getMRVer()->getSSAVersion();

                  auto mrSrcId = graph.addGraphNodeFromMemoryRegion(fun->getName(),mrId,mrVer);

                  graph.addEdge(mrSrcId, instId, "use");
                }
              }
              else if(auto gep = dyn_cast<GepStmt>(stmt)){
                auto src = gep->getRHSVar();
                auto dst = gep->getLHSVar();
                allVarInFun.push_back(dst);
                allVars.push_back(dst);


                auto rhsId = graph.addGraphNodeFromSvfVar(src);
                auto lhsId = graph.addGraphNodeFromSvfVar(dst);
                // todo add edge label.
                graph.addEdge(rhsId, instId, "data");
                graph.addEdge(instId, lhsId, "data");
              }
              else if(auto ret = dyn_cast<RetPE>(stmt)){
                auto src = ret->getRHSVar();
                auto dst = ret->getLHSVar();

                auto rhsId = graph.addGraphNodeFromSvfVar(src);
                auto lhsId = graph.addGraphNodeFromSvfVar(dst);
                // todo add edge label.
                graph.addEdge(rhsId, instId, "data");
                graph.addEdge(instId, lhsId, "data");
              }
              else if(auto copy = dyn_cast<CopyStmt>(stmt)){
                auto src = copy->getRHSVar();
                auto dst = copy->getLHSVar();
                allVarInFun.push_back(dst);
                allVars.push_back(dst);


                auto rhsId = graph.addGraphNodeFromSvfVar(src);
                auto lhsId = graph.addGraphNodeFromSvfVar(dst);
                // todo add edge label.
                graph.addEdge(rhsId, instId, "data");
                graph.addEdge(instId, lhsId, "data");
              }
            }
          }
      }

      for(size_t i = 0; i < allVarInFun.size(); ++i){
        for(size_t j = i+1; j < allVarInFun.size(); ++j){
          auto aliasRes = ander->alias(allVarInFun[i], allVarInFun[j]);
          if(aliasRes == MustAlias){
            // todo: add must alias edge
            auto aliasId1 = graph.addGraphNodeFromSvfVar(allVarInFun[i]);
            auto aliasId2 = graph.addGraphNodeFromSvfVar(allVarInFun[j]);

            graph.addEdge(aliasId1, aliasId2, "mustAlias");
            graph.addEdge(aliasId2, aliasId1, "mustAlias");

          }
          else if(aliasRes == MayAlias || aliasRes == PartialAlias){
            // todo: add may alias edge
            auto aliasId1 = graph.addGraphNodeFromSvfVar(allVarInFun[i]);
            auto aliasId2 = graph.addGraphNodeFromSvfVar(allVarInFun[j]);

            graph.addEdge(aliasId1, aliasId2, "mayAlias");
            graph.addEdge(aliasId2, aliasId1, "mayAlias");
          }
        }
      }

    }

    if(mssa->hasReturnMu(fun)){
      MemSSA::MUSet & return_mus = mssa->getReturnMuSet(fun);
      for(MemSSA::MUSet::iterator mu_it = return_mus.begin(); mu_it != return_mus.end(); mu_it++){
          // no need to process retmu explicitly.
      }
    }
    
  }


  return graph;

}

struct PointsToGraph{
	std::map<u_int32_t, std::map<std::string, std::set<u_int32_t>>> graph;
	u_int32_t nodeNum;
	u_int32_t size;
	std::map<u_int32_t, std::set<u_int32_t>> ptrIdToDerivedMrIdsMap;
	std::map<u_int32_t, std::tuple<std::string, u_int32_t, u_int32_t>> nodeId2mrIdAndVerMap;


	void printRealPtg(){
		// print points-to set
		std::cout << "printing points-to set.\nMax id is: " << nodeNum << "\nTotal entries: " << size << "\n";
		std::vector<std::string> edgeType{"p", "c", "s", "l", "p_d", "p_u", "p_c", "d", "u", "n", "k", "a"};
		for(auto et : edgeType){
			for(u_int32_t i = 0; i < nodeNum; ++i){
				if(graph.count(i) && graph.at(i).count(et)){
					std::cout << i << " == " << et << " ==> {";
					for(auto pointee : graph.at(i).at(et)){
						std::cout << pointee << " ";
					}
					std::cout << "}\n";
				}
			}
		}
	}
	void printNodeId2mrIdAndVerMap(){
		std::cout << "Printing nodeId2mrIdAndVerMap.\n";
		for(auto p : nodeId2mrIdAndVerMap){
			auto nodeId = p.first;
      auto [funName, mrId, mrVer] = p.second;
      std::cout << nodeId << " ==> " << funName << ":MR_" << mrId << "V_" << mrVer << "\n";
    }
	}
	void printPtrIdToDerivedMrIdsMap(){
		std::cout << "Printing ptrIdToDerivedMrIdsMap.\n";
		for(auto p : ptrIdToDerivedMrIdsMap){
        auto nodeId = p.first;
        for(auto mrId : p.second){
            std::cout << "Node id " << nodeId << " => " << "MR id " << mrId << "\n";
        }
    }
	}

};

PointsToGraph readPtgDirectly(const std::set<std::tuple<SVF::NodeID, SVF::NodeID, std::string>> &ptg, 
  const std::map<SVF::NodeID, std::set<SVF::MRID>> &ptrIdToDerivedMrIdsMap, 
  const std::map<SVF::NodeID, std::tuple<std::string, SVF::MRID, SVF::MRVERID>> &nodeId2mrIdAndVerMap){


	PointsToGraph res;
	std::map<u_int32_t, std::map<std::string, std::set<u_int32_t>>> data;
	u_int32_t nodeNum = 0;
	u_int32_t size = 0;

	for(auto [from, to, type] : ptg){
		nodeNum = std::max(nodeNum, from);
		nodeNum = std::max(nodeNum, to);
		data[from][type].insert(to);
		++size;
	}

	res.graph = data;
	res.nodeNum = nodeNum;
	res.size = size;

	res.ptrIdToDerivedMrIdsMap = ptrIdToDerivedMrIdsMap;
	res.nodeId2mrIdAndVerMap = nodeId2mrIdAndVerMap;

	return res;
}




int main(int argc, char **argv) {

  std::vector<std::string> moduleNameVec;
  moduleNameVec = OptionBase::parseOptions(argc, argv, "Whole Program Points-to Analysis", "[options] <input-bitcode...>");

  if(moduleNameVec.size() != 1){
    std::cerr << "Please provide a single llvm-ir file as input.\n";
    exit(1);
  }

  auto graphOutputFileName = getGraphOutputName(moduleNameVec[0]);
  auto outputFileLocation = getGraphOutputLocation(moduleNameVec[0]);

  LLVMModuleSet::buildSVFModule(moduleNameVec);

  /// Build Program Assignment Graph (SVFIR)
  SVFIRBuilder builder;
  SVFIR *pag = builder.build();

  /// Create Andersen's pointer analysis
  Andersen *ander = AndersenWaveDiff::createAndersenWaveDiff(pag);

  SVFGBuilder memSSA;
  auto mssa = memSSA.buildMSSA(ander, true);


  // Initial ptg stored at mssa->ptg.
  mssa->buildInitialPointsToGraph(ander);


  if(Options::DumpIntermediateResultsToFile()){
    dumpNecessaryFilesForGPU(outputFileLocation, graphOutputFileName, mssa.get());
  }


  if(!Options::SbvBased()){
    solve(mssa.get());
    if(Options::DumpFspaPointsToSet()){
      printRealPtg(mssa.get(), false);
    }
  }
  else{
    auto ptgSbv = translatePtgIntoSbv(mssa.get());
    solveSbv(ptgSbv, mssa.get());
    verifyPts(ptgSbv, mssa.get());
    if(Options::DumpFspaPointsToSet()){
      printSbvPts(ptgSbv, mssa.get());
    }
  }
  

  if(Options::DumpIntermediateResultsToFile()){
    dumpGroundTruthPtgToFile(outputFileLocation, graphOutputFileName, mssa.get());
  }


  if(Options::ProGraML()){
    auto graph = generateExtendedProGraMlGraph(mssa.get(), ander, pag);
    std::ofstream graphOut(graphOutputFileName+"-programl.txt");
    if(!graphOut){
      std::cerr << "Cannot open file to dump programl graph\n";
    }
    graphOut << graph.dumpGraph();
    graphOut.close();
  }


  auto ptg = readPtgDirectly(mssa->getTuplePtg(), mssa->getNodeId2MrIdsMap(), mssa->getNodeId2mrIdAndVerMap());
  gpamain(ptg, graphOutputFileName);

  std::cout << "11\n";
  

  return 0;
}
