#ifndef GPA_PROGRAML_H
#define GPA_PROGRAML_H

#include <vector>
#include <string>
#include <set>
#include <cassert>
#include <tuple>
#include <map>
#include "Graphs/ICFG.h"
#include <sstream>
#include "SVFIR/SVFVariables.h"

struct ProGraML{

	// a node is represented as (id, type);
	using NodeType = std::pair<size_t, std::string>;

	public:

		ProGraML(){
			// id for external node is 0.
		}

		size_t addGraphNodeFromIcfgNode(const SVF::ICFGNode *node){
			icfgNodeToGraphIdMap.try_emplace(node, id++);
			idToIcfgNodeMap.try_emplace(icfgNodeToGraphIdMap.at(node), node);
			return icfgNodeToGraphIdMap.at(node);
		}

		size_t addGraphNodeFromSvfVar(const SVF::SVFVar *node){
			svfVarToGraphIdMap.try_emplace(node, id++);
			idToSvfVarMap.try_emplace(svfVarToGraphIdMap.at(node), node);
			return svfVarToGraphIdMap.at(node);
		}

		size_t addGraphNodeFromMemoryRegion(std::string funName, size_t mrId, size_t verId){
			std::tuple<std::string, size_t, size_t> t(funName, mrId, verId);
			memoryRegionToGraphIdMap.try_emplace(t, id++);
			idToMemoryRegionMap.try_emplace(memoryRegionToGraphIdMap.at(t), t);
			return memoryRegionToGraphIdMap.at(t);
		}

		void addEdge(size_t from, size_t to, std::string edgeType){
			assert(nodesToInfoMap.count(from) && "Unknown from node.");
			assert(nodesToInfoMap.count(to) && "Unknown to node.");

			edges[from][to].insert(edgeType);

		}

		std::string dumpGraph(){
			std::string str;
      std::stringstream ss(str); 

			ss << "idToIcfgNodeMap\n";
			for(auto p : idToIcfgNodeMap){
				ss << p.first << " => ICFGNode: " << p.second->getId() << "\n";
			}

			ss << "idToSvfVarMap\n";
			for(auto p : idToSvfVarMap){
				ss << p.first << " => SVFVar: " << p.second->getId() << "\n";
			}

			ss << "idToMemoryRegionMap\n";
			for(auto p : idToMemoryRegionMap){
				auto mr = p.second;
				ss << p.first << " => MR: " << std::get<0>(mr) << "_" + std::to_string(std::get<1>(mr)) << "_" + std::to_string(std::get<2>(mr)) << "\n";
			}

			for(auto p0 : edges){
				auto from = p0.first;
				for(auto p1 : p0.second){
					auto to = p1.first;
					for(auto e : p1.second){
						ss << from << " == " << e << " ==> " << to << "\n";
					}
				}
			}

			for(size_t i = 0; i < id; ++i){
				if(i == 0){
					ss << "G.node(\"" << i << "\", label=\"external\")\n";
				}
				else if(!getId(i).empty()){
					ss << "G.node(\"" << i << "\", label=\"" << getId(i) << "\")\n";
				}
				
			}

			std::map<std::string, std::string> colors{{"call", "green"}, {"data", "red"}, {"control", "black"},
																								{"def", "blue"}, {"use", "yellow"}, {"mustAlias", "cyan"},
																								{"mayAlias", "purple"}, {"test", "orange"}};

			for(auto p0 : edges){
				auto from = p0.first;
				for(auto p1 : p0.second){
					auto to = p1.first;
					for(auto e : p1.second){
						ss << "G.edge(\"" << from << "\", \"" << to << "\", label=\"" << e << "\", color=\"" 
								<< colors[e] << "\", penwidth='3.0')\n";
					}
				}
			}
      return ss.str();
		}




	private:

		inline static size_t id = 1;

		std::map<size_t, std::string> nodesToInfoMap;
		std::map<size_t, std::map<size_t, std::set<std::string>>> edges;
		std::map<const SVF::ICFGNode*, size_t> icfgNodeToGraphIdMap;
		std::map<const SVF::SVFVar*, size_t> svfVarToGraphIdMap;
		std::map<std::tuple<std::string, size_t, size_t>, size_t> memoryRegionToGraphIdMap;
		std::map<size_t, const SVF::ICFGNode*> idToIcfgNodeMap;
		std::map<size_t, const SVF::SVFVar*> idToSvfVarMap;
		std::map<size_t, std::tuple<std::string, size_t, size_t>> idToMemoryRegionMap;

		std::string getId(size_t nodeId){
			if(idToIcfgNodeMap.count(nodeId)){
				return "ICFGNode: " + std::to_string(idToIcfgNodeMap.at(nodeId)->getId());
			}
			else if(idToSvfVarMap.count(nodeId)){
				return "SVFVar: " + std::to_string(idToSvfVarMap.at(nodeId)->getId());
			}
			else if(idToMemoryRegionMap.count(nodeId)){
				auto mr = idToMemoryRegionMap.at(nodeId);
				return "MR: " + std::get<0>(mr) + "_" + std::to_string(std::get<1>(mr)) + "_" + std::to_string(std::get<2>(mr));
			}
			else{
				return "";
			}
		}




};


#endif