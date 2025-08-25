#include "PointsToGraph.h"

#include "MSSA/MemSSA.h"


PointsToGraph::PointsToGraph(SVF::MemSSA *mssa){
	std::map<u_int32_t, std::map<std::string, std::set<u_int32_t>>> data;
	u_int32_t nodeNum = 0;
	u_int32_t size = 0;

	for(auto [from, to, type] : mssa->getTuplePtg()){
		nodeNum = std::max(nodeNum, from);
		nodeNum = std::max(nodeNum, to);
		data[from][type].insert(to);
		++size;
	}

	this->graph = data;
	this->nodeNum = nodeNum;
	this->size = size;

	this->ptrIdToDerivedMrIdsMap = mssa->getNodeId2MrIdsMap();
	this->nodeId2mrIdAndVerMap = mssa->getNodeId2mrIdAndVerMap();
}


void PointsToGraph::printRealPtg(){
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

void PointsToGraph::printNodeId2mrIdAndVerMap(){
	std::cout << "Printing nodeId2mrIdAndVerMap.\n";
	for(auto p : nodeId2mrIdAndVerMap){
		auto nodeId = p.first;
		auto [funName, mrId, mrVer] = p.second;
		std::cout << nodeId << " ==> " << funName << ":MR_" << mrId << "V_" << mrVer << "\n";
	}
}
void PointsToGraph::printPtrIdToDerivedMrIdsMap(){
	std::cout << "Printing ptrIdToDerivedMrIdsMap.\n";
	for(auto p : ptrIdToDerivedMrIdsMap){
			auto nodeId = p.first;
			for(auto mrId : p.second){
					std::cout << "Node id " << nodeId << " => " << "MR id " << mrId << "\n";
			}
	}
}