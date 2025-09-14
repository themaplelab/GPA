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

		ProGraML(SVF::SVFIR *pag, std::map<u_int32_t, std::set<u_int32_t>> ptrIdToDerivedMrIdsMap) : pag(pag){
			// id for external node is 0.
			for(auto pair : ptrIdToDerivedMrIdsMap){
				auto svfNodeId = pair.first;
				for(auto mrId : pair.second){
					mrIdToSvfNodeIdMap[mrId].insert(svfNodeId);
				}
			}

			
		}

		size_t addGraphNodeFromIcfgNode(const SVF::ICFGNode *node){

			auto [it, inserted] = icfgNodeToGraphIdMap.try_emplace(node, id);
			if(inserted){
				idToIcfgNodeMap.emplace(id, node);
				++id;
			}
			return it->second;


			// icfgNodeToGraphIdMap.try_emplace(node, id++);
			// idToIcfgNodeMap.try_emplace(icfgNodeToGraphIdMap.at(node), node);
			// return icfgNodeToGraphIdMap.at(node);
		}

		size_t addGraphNodeFromSvfVar(const SVF::SVFVar *node){

			auto [it, inserted] = svfVarToGraphIdMap.try_emplace(node, id);
			if(inserted){
				idToSvfVarMap.emplace(id, node);
				++id;
			}
			return it->second;


			// svfVarToGraphIdMap.try_emplace(node, id++);
			// idToSvfVarMap.try_emplace(svfVarToGraphIdMap.at(node), node);
			// return svfVarToGraphIdMap.at(node);
		}

		size_t addGraphNodeFromMemoryRegion(std::string funName, size_t mrId, size_t verId){


			std::tuple<std::string, size_t, size_t> t(funName, mrId, verId);
			auto [it, inserted] = memoryRegionToGraphIdMap.try_emplace(t, id);
			if(inserted){
				idToMemoryRegionMap.emplace(id, t);
				++id;
			}
			return it->second;

			// std::tuple<std::string, size_t, size_t> t(funName, mrId, verId);
			// memoryRegionToGraphIdMap.try_emplace(t, id++);
			// idToMemoryRegionMap.try_emplace(memoryRegionToGraphIdMap.at(t), t);
			// return memoryRegionToGraphIdMap.at(t);
		}

		void addEdge(size_t from, size_t to, std::string edgeType, size_t pos = -1){
			assert(nodesToInfoMap.count(from) && "Unknown from node.");
			assert(nodesToInfoMap.count(to) && "Unknown to node.");

			edges[from][to].insert(edgeType);
			if(edgeType == "control" || edgeType == "data"){
				std::tuple<size_t, size_t, std::string> t{from, to, edgeType};
				positionsMap.try_emplace(t, pos);
			}

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

			for(auto pair : positionsMap){
				auto [from, to, et] = pair.first;
				ss << "(" << from << ", " << et << ", " << to << ") pos: " << pair.second << "\n";
			}

			for(size_t i = 1; i < id; ++i){
				

				if(idToIcfgNodeMap.count(i)){

				}
				else if(idToSvfVarMap.count(i)){
					ss << "ProGraMl node id: " << i << " has type " << idToSvfVarMap.at(i)->getType()->getKind() << "\n";
				}
				else if(idToMemoryRegionMap.count(i)){
					auto svfNodeIds = mrIdToSvfNodeIdMap.at(std::get<1>(idToMemoryRegionMap.at(i)));
					for(auto sni : svfNodeIds){
						ss << "ProGraMl node id: " << i << " has type " << pag->getGNode(sni)->getType()->getKind() << "\n";
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

		std::string dumpDataSet(const std::string &benchmarkName){
			std::string str;
      std::stringstream ss(str); 

			// dump metadata
			ss << benchmarkName << "," << id << "," << edges.size() << "\n";

			// dump nodes
			
			for(size_t i = 0; i < id; ++i){
				if(i == 0){
					ss << i << "," << "Instruction" << "," << "NA" << "\n";
				}
				else{
					ss << i << "," << getKind(i) << "," << getDtype(i) << "\n";
				}
			}


			// dump edges

			for(auto p0 : edges){
				auto from = p0.first;
				for(auto p1 : p0.second){
					auto to = p1.first;
					for(auto e : p1.second){
						std::tuple<size_t, size_t, std::string> t{from, to, e};
						ss << from << "," << to << "," << e << "," << getPos(t) << "\n";
					}
				}
			}
		}


	private:

		inline static size_t id = 1;

		std::map<size_t, std::string> nodesToInfoMap;
		std::map<size_t, std::map<size_t, std::set<std::string>>> edges;
		std::map<std::tuple<size_t, size_t, std::string>, size_t> positionsMap;
		std::map<const SVF::ICFGNode*, size_t> icfgNodeToGraphIdMap;
		std::map<const SVF::SVFVar*, size_t> svfVarToGraphIdMap;
		std::map<std::tuple<std::string, size_t, size_t>, size_t> memoryRegionToGraphIdMap;
		std::map<size_t, const SVF::ICFGNode*> idToIcfgNodeMap;
		std::map<size_t, const SVF::SVFVar*> idToSvfVarMap;
		std::map<size_t, std::tuple<std::string, size_t, size_t>> idToMemoryRegionMap;
		std::map<u_int32_t, std::set<u_int32_t>> mrIdToSvfNodeIdMap;;
		// std::map<const SVF::SVFVar*, SVF::SVFType::SVFTyKind> svfVarToTypeMap;

		SVF::SVFIR *pag;

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

		std::string getKind(size_t nodeId){
			if(idToIcfgNodeMap.count(nodeId)){
				return "Instruction";
			}
			else if(idToSvfVarMap.count(nodeId)){
				return "Top-level variable";
			}
			else if(idToMemoryRegionMap.count(nodeId)){
				auto mr = idToMemoryRegionMap.at(nodeId);
				return "Address-taken variable";
			}
			else{
				return "";
			}
		}

		std::string getDtype(size_t nodeId){
			if(idToIcfgNodeMap.count(nodeId)){
				return "";
			}
			else if(idToSvfVarMap.count(nodeId)){
				return std::to_string(idToSvfVarMap.at(nodeId)->getType()->getKind());
			}
			else if(idToMemoryRegionMap.count(nodeId)){
				auto svfNodeIds = mrIdToSvfNodeIdMap.at(std::get<1>(idToMemoryRegionMap.at(nodeId)));
				for(auto sni : svfNodeIds){
					return std::to_string(pag->getGNode(sni)->getType()->getKind());
				}
			}
			else{
				return "";
			}
		}

		std::string getPos(std::tuple<size_t, size_t, std::string> t){
			if(positionsMap.count(t)){
				return std::to_string(positionsMap.at(t));
			}
			else{
				return "NA";
			}
		}


};


#endif