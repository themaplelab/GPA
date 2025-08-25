#ifndef GPA_POINTS_TO_GRAPH_H
#define GPA_POINTS_TO_GRAPH_H


#include <map>
#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <string>


namespace SVF {
    class MemSSA;   // forward declaration
}



struct PointsToGraph{
	std::map<u_int32_t, std::map<std::string, std::set<u_int32_t>>> graph;
	u_int32_t nodeNum;
	u_int32_t size;
	std::map<u_int32_t, std::set<u_int32_t>> ptrIdToDerivedMrIdsMap;
	std::map<u_int32_t, std::tuple<std::string, u_int32_t, u_int32_t>> nodeId2mrIdAndVerMap;


	explicit PointsToGraph(SVF::MemSSA *mssa);
	PointsToGraph() = default;

	void printRealPtg();
	void printNodeId2mrIdAndVerMap();
	void printPtrIdToDerivedMrIdsMap();

};



#endif