#include "gpa.cuh"
#include <cstdio>


#include <iostream>
#include <tuple>
#include <fstream>
#include <set>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cassert>
#include <chrono>

#include "sbv_gpu.cuh"
// #include "test.cuh"


#define DEBUG_LOG(x) do { if (debug) std::cout << x << std::endl; } while(0)

#define MAX_EDGE_BUFFER_SIZE 4096



// int test(){
	
// 	sbvInitUnitTest();
// 	sbvSetUnitTest();
// 	sbvContainUnitTest();
// 	sbvToArrayUnitTest();

// 	return 0;

// }

struct PointerSolver{


};

struct PointsToGraph{
	std::map<u_int32_t, std::map<std::string, std::set<u_int32_t>>> graph;
	u_int32_t nodeNum;
	u_int32_t size;
	std::map<u_int32_t, std::set<u_int32_t>> ptrIdToDerivedMrIdsMap;
	std::map<u_int32_t, std::tuple<std::string, u_int32_t, u_int32_t>> nodeId2mrIdAndVerMap;
	PointerSolver solver;


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

PointsToGraph readPtg(std::string fileName){
	std::ifstream inFile(fileName+"-initialPtg.txt");
	if(!inFile){
		std::cerr << "Cannot open file " << fileName << "-initialPtg.txt\n";
		std::terminate();
	}

	PointsToGraph res;
	std::map<u_int32_t, std::map<std::string, std::set<u_int32_t>>> data;
	u_int32_t nodeNum = 0;
	u_int32_t size = 0;


	std::string line;
	while(std::getline(inFile, line)){
		std::stringstream ss(line);
		std::string element;
		std::vector<std::string> elements;
		while(std::getline(ss, element, ',')){
			elements.push_back(element);
		}
		u_int32_t from = static_cast<uint32_t>(std::stoul(elements[0]));
		u_int32_t to = static_cast<uint32_t>(std::stoul(elements[1]));
		nodeNum = std::max(nodeNum, from);
		nodeNum = std::max(nodeNum, to);
		data[from][elements[2]].insert(to);
		++size;
	}

	res.graph = data;
	res.nodeNum = nodeNum;
	res.size = size;

	inFile.close();


	std::map<u_int32_t, std::set<u_int32_t>> ptrIdToDerivedMrIdsMap;
	std::ifstream ptrIdToDerivedMrIdsMapInFile(fileName+"-nodeId2RelatedMrId.txt");
	if(!ptrIdToDerivedMrIdsMapInFile){
		std::cerr << "Cannot open file " << fileName << "-nodeId2RelatedMrId.txt\n";
		std::terminate();
	}

	while(std::getline(ptrIdToDerivedMrIdsMapInFile, line)){
		std::stringstream ss(line);
		std::string element;
		std::vector<std::string> elements;
		while(std::getline(ss, element, ',')){
			elements.push_back(element);
		}
		u_int32_t ptrId = static_cast<uint32_t>(std::stoul(elements[0]));
		std::set<u_int32_t> mrIds;
		for(int i = 1; i < elements.size(); ++i){
			mrIds.insert(static_cast<uint32_t>(std::stoul(elements[i])));
		}
		
		ptrIdToDerivedMrIdsMap.emplace(ptrId, mrIds);
	}

	res.ptrIdToDerivedMrIdsMap = ptrIdToDerivedMrIdsMap;
	ptrIdToDerivedMrIdsMapInFile.close();


	std::map<u_int32_t, std::tuple<std::string, u_int32_t, u_int32_t>> nodeId2mrIdAndVerMap;
	std::ifstream nodeId2mrIdAndVerMapInFile(fileName+"-nodeId2MrId.txt");
	if(!nodeId2mrIdAndVerMapInFile){
		std::cerr << "Cannot open file " << fileName << "-nodeId2MrId.txt\n";
		std::terminate();
	}

	while(std::getline(nodeId2mrIdAndVerMapInFile, line)){
		std::stringstream ss(line);
		std::string element;
		std::vector<std::string> elements;
		while(std::getline(ss, element, ',')){
			elements.push_back(element);
		}

		u_int32_t nodeId = static_cast<uint32_t>(std::stoul(elements[0]));
		uint32_t mrId = static_cast<uint32_t>(std::stoul(elements[2]));
		uint32_t verId = static_cast<uint32_t>(std::stoul(elements[3]));

		std::tuple<std::string, u_int32_t, u_int32_t> t(elements[1], mrId, verId);
		nodeId2mrIdAndVerMap.emplace(nodeId, t);
	}

	res.nodeId2mrIdAndVerMap = nodeId2mrIdAndVerMap;
	nodeId2mrIdAndVerMapInFile.close();

	return res;
}


std::map<std::string, std::vector<SBVCPU*>> preprocessPtgForGpu(PointsToGraph &ptg){
	auto graph = ptg.graph;
	std::vector<std::string> edgeType{"p", "c", "s", "l", "p_d", "p_u", "p_c", "d", "u", "n", "k", "a"};

	std::map<std::string, std::vector<SBVCPU*>> totalSbvs;
	for(auto et : edgeType){
		std::vector<SBVCPU*> sbvs;
		// 0-th entry is intentially introduced because the nodeId starts from 1.
		for(int i = 0; i <= ptg.nodeNum; ++i){
			sbvs.push_back(new SBVCPU());
		}
		totalSbvs.emplace(et, sbvs);
	}

	for(auto pair : graph){
		auto from = pair.first;
		for(auto typeAndTos : pair.second){
			auto type = typeAndTos.first;
			for(auto to : typeAndTos.second){
				totalSbvs[type][from]->setCPU(to);
			}
		}
	}

	return totalSbvs;
}

std::map<std::string, std::vector<SBVCPU*>> preprocessPtgForGpuWithReverseEdges(PointsToGraph &ptg){
	auto graph = ptg.graph;
	std::set<std::string> forwardEdgeType{"p", "u", "s", "p_c", "p_u", "u", "a", "k"};
	std::set<std::string> reverseEdgeType{"c", "l", "d", "p_d", "k", "n"};


	std::map<std::string, std::vector<SBVCPU*>> totalSbvs;
	for(auto et : forwardEdgeType){
		std::vector<SBVCPU*> sbvs;
		// 0-th entry is intentially introduced because the nodeId starts from 1.
		for(int i = 0; i <= ptg.nodeNum; ++i){
			sbvs.push_back(new SBVCPU());
		}
		totalSbvs.emplace(et, sbvs);
	}

	for(auto et : reverseEdgeType){
		std::vector<SBVCPU*> sbvs;
		// 0-th entry is intentially introduced because the nodeId starts from 1.
		for(int i = 0; i <= ptg.nodeNum; ++i){
			sbvs.push_back(new SBVCPU());
		}
		totalSbvs.emplace(et+"-1", sbvs);
	}

	for(auto pair : graph){
		auto from = pair.first;
		for(auto typeAndTos : pair.second){
			auto type = typeAndTos.first;

			if(type == "k"){
				for(auto to : typeAndTos.second){
					totalSbvs[type][from]->setCPU(to);
					totalSbvs[type+"-1"][to]->setCPU(from);
				}
			}
			else if(forwardEdgeType.count(type)){
				for(auto to : typeAndTos.second){
					totalSbvs[type][from]->setCPU(to);
				}
			}
			else{
				for(auto to : typeAndTos.second){
					totalSbvs[type+"-1"][to]->setCPU(from);
				}
			}
			
		}
	}

	// for(auto p : totalSbvs){
	// 	std::cout << p.first << "\n";
	// }

	return totalSbvs;
}



SBV* gpuMallocSbvFromHost(SBV* pool, u_int32_t* pool_ptr) {
    SBV* result = pool + (*pool_ptr);  // just pointer math, no dereferencing
    (*pool_ptr)++;
    return result;
}

__device__ u_int32_t getMrIdOfPtgNode(u_int32_t node, const u_int32_t* nodeIdToMrId, u_int32_t maxNodeId){
		if(node > maxNodeId){
			return 0;
		} 
    return nodeIdToMrId[node];
}

__device__ bool isaVersionofPtr(u_int32_t ptrId, u_int32_t mrAndVerId, const u_int32_t* nodeIdToMrId, u_int32_t maxNodeId,
    const u_int32_t* derived_ids, const u_int32_t* derived_offsets, const u_int32_t* derived_counts){
	
			// ptrid 17 meAndVerId 65

	u_int32_t mrId = getMrIdOfPtgNode(mrAndVerId, nodeIdToMrId, maxNodeId);
	u_int32_t offset = derived_offsets[ptrId];
	u_int32_t count = derived_counts[ptrId];

	// printf("%d %d %d %d %d\n", mrId, offset, count, ptrId, mrAndVerId);

	for (u_int32_t i = 0; i < count; ++i) {
			if (derived_ids[offset + i] == mrId) {
					return true;
			}
	}
	return false;
}

__global__ void rule1(SBV_ADDR_TYPE *pData, SBV_ADDR_TYPE *cData, u_int32_t nodeNum, SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed){
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// valid id range from [1, nodeNum]
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE cAddr = cData[idx];
	SBV_ADDR_TYPE pAddr = pData[idx];
	if(cAddr == NULL_ADDR || pAddr == NULL_ADDR){
		return;
	}

	SBV *cDstSbv = &sbvPool[cAddr];
	SBV *pDstSbv = &sbvPool[pAddr];
	if(cDstSbv->empty(sbvPool) || pDstSbv->empty(sbvPool)){
		return;
	}

	u_int32_t cbuffer[MAX_EDGE_BUFFER_SIZE];
	int clen = cDstSbv->toArray(cbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	u_int32_t pbuffer[MAX_EDGE_BUFFER_SIZE];
	int plen = pDstSbv->toArray(pbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < clen; ++i){
		auto cDst = cbuffer[i];
		SBV_ADDR_TYPE pUpdateAddr = pData[cDst];
		for(int j = 0; j < plen; ++j){
			auto pDst = pbuffer[j];
			// printf("trying Add p-edge (%d,%d) %d\n", cDst, pDst, idx);
			if(sbvPool[pUpdateAddr].set(pDst, pUpdateAddr, sbvPool, sbvPoolIndex)){

				printf("Adding p-edge (%d,%d)\n", cDst, pDst);
				atomicExch(changed, 1);
			}
		}
	}
}

__global__ void rule2(SBV_ADDR_TYPE *uData, SBV_ADDR_TYPE *lData, SBV_ADDR_TYPE *pcData, SBV_ADDR_TYPE *cData, u_int32_t nodeNum, 
										SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE uAddr = uData[idx];
	SBV_ADDR_TYPE lAddr = lData[idx];

	if(uAddr == NULL_ADDR || lAddr == NULL_ADDR){
		return;
	}


	SBV *uDstSbv = &sbvPool[uAddr];
	SBV *lDstSbv = &sbvPool[lAddr];
	if(uDstSbv->empty(sbvPool) || lDstSbv->empty(sbvPool)){
		return;
	}

	u_int32_t ubuffer[MAX_EDGE_BUFFER_SIZE];
	int ulen = uDstSbv->toArray(ubuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	u_int32_t lbuffer[MAX_EDGE_BUFFER_SIZE];
	int llen = lDstSbv->toArray(lbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < ulen; ++i){
		auto uDst = ubuffer[i];
		SBV_ADDR_TYPE pcDstAddr = pcData[uDst];
		if(sbvPool[pcDstAddr].empty(sbvPool)){
			continue;
		}
		SBV_ADDR_TYPE cUpdateAddr = cData[uDst];
		for(int j = 0; j < llen; ++j){
			auto lDst = lbuffer[j];
			if(sbvPool[pcDstAddr].test(lDst, pcDstAddr, sbvPool)){
				if(sbvPool[cUpdateAddr].set(lDst, cUpdateAddr, sbvPool, sbvPoolIndex)){
					printf("Adding c-edge (%d,%d)\n", uDst, lDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
	return;
}

__global__ void rule3(SBV_ADDR_TYPE *dData, SBV_ADDR_TYPE *sData, SBV_ADDR_TYPE *pcData, SBV_ADDR_TYPE *cData, u_int32_t nodeNum, 
										SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE dAddr = dData[idx];
	SBV_ADDR_TYPE sAddr = sData[idx];

	if(dAddr == NULL_ADDR || sAddr == NULL_ADDR){
		return;
	}

	SBV *dDstSbv = &sbvPool[dAddr];
	SBV *sDstSbv = &sbvPool[sAddr];
	if(dDstSbv->empty(sbvPool) || sDstSbv->empty(sbvPool)){
		return;
	}

	u_int32_t sbuffer[MAX_EDGE_BUFFER_SIZE];
	int slen = sDstSbv->toArray(sbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	u_int32_t dbuffer[MAX_EDGE_BUFFER_SIZE];
	int dlen = dDstSbv->toArray(dbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < slen; ++i){
		auto sDst = sbuffer[i];
		 
		SBV_ADDR_TYPE pcDstAddr = pcData[sDst];
		if(sbvPool[pcDstAddr].empty(sbvPool)){
			continue;
		}

		
		SBV_ADDR_TYPE cUpdateAddr = cData[sDst];
		for(int j = 0; j < dlen; ++j){
			auto dDst = dbuffer[j];
			if(sbvPool[pcDstAddr].test(dDst, pcDstAddr, sbvPool)){
				if(sbvPool[cUpdateAddr].set(dDst, cUpdateAddr, sbvPool, sbvPoolIndex)){
					printf("Adding c-edge (%d,%d)\n", sDst, dDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
}

__global__ void rule4(SBV_ADDR_TYPE *pdData, SBV_ADDR_TYPE *pData, SBV_ADDR_TYPE *dData, u_int32_t nodeNum, SBV *sbvPool, 
					u_int32_t *sbvPoolIndex, int *changed, const u_int32_t* nodeIdToMrId,
    			const u_int32_t* derivedIds, const u_int32_t* derivedOffsets, const u_int32_t* derivedCounts){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE pdAddr = pdData[idx];
	SBV_ADDR_TYPE pAddr = pData[idx];
	if(pdAddr == NULL_ADDR || pAddr == NULL_ADDR){
		return;
	}

	SBV *pdDstSbv = &sbvPool[pdAddr];
	SBV *pDstSbv = &sbvPool[pAddr];
	if(pdDstSbv->empty(sbvPool) || pDstSbv->empty(sbvPool)){
		return;
	}

	SBV_ADDR_TYPE dDstAddr = dData[idx];
	u_int32_t pdbuffer[MAX_EDGE_BUFFER_SIZE];
	int pdlen = pdDstSbv->toArray(pdbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	u_int32_t pbuffer[MAX_EDGE_BUFFER_SIZE];
	int plen = pDstSbv->toArray(pbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < pdlen; ++i){
		auto pdDst = pdbuffer[i];
		for(int j = 0; j < plen; ++j){
			auto pDst = pbuffer[j];
			if(isaVersionofPtr(pDst, pdDst, nodeIdToMrId, nodeNum, derivedIds, derivedOffsets, derivedCounts)){
				if(sbvPool[dDstAddr].set(pdDst, dDstAddr, sbvPool, sbvPoolIndex)){
					printf("Adding d-edge (%d,%d)\n", idx, pdDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
	return;
}

__global__ void rule5(SBV_ADDR_TYPE *puData, SBV_ADDR_TYPE *pData, SBV_ADDR_TYPE *uData, u_int32_t nodeNum, SBV *sbvPool, 
					u_int32_t *sbvPoolIndex, int *changed, const u_int32_t* nodeIdToMrId,
    			const u_int32_t* derivedIds, const u_int32_t* derivedOffsets, const u_int32_t* derivedCounts){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE puDstAddr = puData[idx];
	SBV_ADDR_TYPE pDstAddr = pData[idx];

	if(puDstAddr == NULL_ADDR || pDstAddr == NULL_ADDR){
		return;
	}


	SBV *puDstSbv = &sbvPool[puDstAddr];
	SBV *pDstSbv = &sbvPool[pDstAddr];
	if(puDstSbv->empty(sbvPool) || pDstSbv->empty(sbvPool)){
		return;
	}

	SBV_ADDR_TYPE uDstAddr = uData[idx];
	u_int32_t pubuffer[MAX_EDGE_BUFFER_SIZE];
	int pulen = puDstSbv->toArray(pubuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	u_int32_t pbuffer[MAX_EDGE_BUFFER_SIZE];
	int plen = pDstSbv->toArray(pbuffer,sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < pulen; ++i){
		auto puDst = pubuffer[i];
		for(int j = 0; j < plen; ++j){
			auto pDst = pbuffer[j];
			// printf("(%d,%d,%d)\n", pDst, puDst, isaVersionofPtr(pDst, puDst, nodeIdToMrId, nodeNum, derivedIds, derivedOffsets, derivedCounts));
			if(isaVersionofPtr(pDst, puDst, nodeIdToMrId, nodeNum, derivedIds, derivedOffsets, derivedCounts)){
				if(sbvPool[uDstAddr].set(puDst, uDstAddr, sbvPool, sbvPoolIndex)){
					printf("Adding u-edge (%d,%d)\n", idx, puDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
}

__global__ void rule6(SBV_ADDR_TYPE *dData, SBV_ADDR_TYPE *kData, SBV_ADDR_TYPE *aData, u_int32_t nodeNum, SBV *sbvPool, 
					u_int32_t *sbvPoolIndex, int *changed){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}



	SBV_ADDR_TYPE dDstAddr = dData[idx];
	SBV *dDstSbv = &sbvPool[dDstAddr];
	if(dDstSbv->empty(sbvPool)){
		return;
	}

	u_int32_t dbuffer[MAX_EDGE_BUFFER_SIZE];
	int dlen = dDstSbv->toArray(dbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	for(int i = 0; i < dlen; ++i){
		auto dDst = dbuffer[i];
		
		SBV_ADDR_TYPE kDstAddr = kData[dDst];
		SBV *kDstSbv = &sbvPool[kDstAddr];
		if(kDstSbv->empty(sbvPool)){
			continue;
		}

		u_int32_t kbuffer[MAX_EDGE_BUFFER_SIZE];
		int klen = kDstSbv->toArray(kbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
		for(int j = 0; j < klen; ++j){
			auto kDst = kbuffer[j];
			
			SBV_ADDR_TYPE aDstAddr = aData[kDst];
			SBV *aDstSbv = &sbvPool[aDstAddr];
			if(aDstSbv->set(dDst, aDstAddr, sbvPool, sbvPoolIndex)){
				printf("Adding a-edge (%d,%d)\n", kDst, dDst);
				atomicExch(changed, 1);
			}
		}
	}
}

__global__ void rule7(SBV_ADDR_TYPE *nData, SBV_ADDR_TYPE *kData, SBV_ADDR_TYPE *aData, SBV_ADDR_TYPE *cData, u_int32_t nodeNum, 
					SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed, const u_int32_t* nodeIdToMrId){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE nDstAddr = nData[idx];
	if(nDstAddr == NULL_ADDR){
		return;
	}
	SBV *nDstSbv = &sbvPool[nDstAddr];
	if(nDstSbv->empty(sbvPool)){
		return;
	}


	SBV_ADDR_TYPE cDstAddr = cData[idx];
	const u_int32_t selfMrId = getMrIdOfPtgNode(idx, nodeIdToMrId, nodeNum);

	u_int32_t nbuffer[MAX_EDGE_BUFFER_SIZE];
	int nlen = nDstSbv->toArray(nbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	for(int i = 0; i < nlen; ++i){
		auto nDst = nbuffer[i];
		
		SBV_ADDR_TYPE kDstAddr = kData[nDst];
		SBV *kDstSbv = &sbvPool[kDstAddr];
		if(kDstSbv->empty(sbvPool)){
			continue;
		}
		u_int32_t kbuffer[MAX_EDGE_BUFFER_SIZE];
		int klen = kDstSbv->toArray(kbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
		for(int j = 0; j < klen; ++j){
			auto kDst = kbuffer[j];
			SBV_ADDR_TYPE aDstAddr = aData[kDst];
			SBV *aDstSbv = &sbvPool[aDstAddr];
			if(aDstSbv->empty(sbvPool)){
				continue;
			}
			u_int32_t abuffer[MAX_EDGE_BUFFER_SIZE];
			int alen = aDstSbv->toArray(abuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
			for(int k = 0; k < alen; ++k){
				auto aDst = abuffer[k];
				if(selfMrId != getMrIdOfPtgNode(aDst, nodeIdToMrId, nodeNum)){
					if(sbvPool[cDstAddr].set(nDst, cDstAddr, sbvPool, sbvPoolIndex)){
						printf("Adding c-edge (%d,%d)\n", idx, nDst);
						atomicExch(changed, 1);
					}
				}
			}
		}
	}
}

void solve(std::map<std::string, SBV_ADDR_TYPE*> dataMap, u_int32_t nodeNum, SBV *deviceSbvPool, u_int32_t *deviceSbvPoolIndex,
						const u_int32_t* nodeIdToMrId, const u_int32_t* derivedIds, const u_int32_t* derivedOffsets, 
						const u_int32_t* derivedCounts, bool debug = false, size_t threadsPerBlock = 16){

	int hostChanged = 1;
	int *deviceChanged;
	cudaMalloc(&deviceChanged, sizeof(int));

	size_t round = 0;
	// valid node id start from 1.
	size_t blockPerGrid = (nodeNum+threadsPerBlock) / threadsPerBlock;
	std::cout << blockPerGrid << " " << threadsPerBlock << "\n";
	while(hostChanged){
		++round;
		
		hostChanged = 0;
		cudaMemcpy(deviceChanged, &hostChanged, sizeof(int), cudaMemcpyHostToDevice);

		rule1<<<blockPerGrid,threadsPerBlock>>>(dataMap["p"], dataMap["c"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		CUDA_CHECK(cudaGetLastError());
		rule2<<<blockPerGrid,threadsPerBlock>>>(dataMap["u"], dataMap["l"], dataMap["p_c"], dataMap["c"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		CUDA_CHECK(cudaGetLastError());
		rule3<<<blockPerGrid,threadsPerBlock>>>(dataMap["d"], dataMap["s"], dataMap["p_c"], dataMap["c"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		CUDA_CHECK(cudaGetLastError());
		rule4<<<blockPerGrid,threadsPerBlock>>>(dataMap["p_d"], dataMap["p"], dataMap["d"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged, nodeIdToMrId, derivedIds,
										derivedOffsets, derivedCounts);
		CUDA_CHECK(cudaGetLastError());
		rule5<<<blockPerGrid,threadsPerBlock>>>(dataMap["p_u"], dataMap["p"], dataMap["u"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged, nodeIdToMrId, derivedIds,
										derivedOffsets, derivedCounts);
		CUDA_CHECK(cudaGetLastError());
		rule6<<<blockPerGrid,threadsPerBlock>>>(dataMap["d"], dataMap["k"], dataMap["a"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		CUDA_CHECK(cudaGetLastError());
		rule7<<<blockPerGrid,threadsPerBlock>>>(dataMap["n"], dataMap["k"], dataMap["a"], dataMap["c"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged, nodeIdToMrId);

		CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(&hostChanged, deviceChanged, sizeof(int), cudaMemcpyDeviceToHost));
	}
	DEBUG_LOG("Round: " << round);

}

__global__ void rule1Reverse(SBV_ADDR_TYPE *pData, SBV_ADDR_TYPE *cDataReverse, u_int32_t nodeNum, SBV *sbvPool, 
					u_int32_t *sbvPoolIndex, int *changed){
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// valid id range from [1, nodeNum]
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE cDstAddr = cDataReverse[idx];
	if(cDstAddr == NULL_ADDR){
		return;
	}

	SBV *cDstSbv = &sbvPool[cDstAddr];
	if(cDstSbv->empty(sbvPool)){
		return;
	}

	u_int32_t cbuffer[MAX_EDGE_BUFFER_SIZE];
	int clen = cDstSbv->toArray(cbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	SBV_ADDR_TYPE pUpdateAddr = pData[idx];
	if(pUpdateAddr == NULL_ADDR){
		return;
	}
	SBV *pUpdateSbv = &sbvPool[pUpdateAddr];

	for(int i = 0; i < clen; ++i){
		auto cDst = cbuffer[i];

		SBV_ADDR_TYPE pDstAddr = pData[cDst];
		if(pDstAddr == NULL_ADDR){
			return;
		}
		SBV *pDstSbv = &sbvPool[pDstAddr];
		u_int32_t pbuffer[MAX_EDGE_BUFFER_SIZE];
		int plen = pDstSbv->toArray(pbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
		for(int j = 0; j < plen; ++j){
			auto pDst = pbuffer[j];
			// printf("trying Add p-edge (%d,%d) %d\n", cDst, pDst, idx);
			if(pUpdateSbv->set(pDst, pUpdateAddr, sbvPool, sbvPoolIndex)){
				// printf("Adding p-edge (%d,%d)\n", idx, pDst);
				atomicExch(changed, 1);
			}
		}
	}
}

__global__ void rule2Reverse(SBV_ADDR_TYPE *uData, SBV_ADDR_TYPE *lDataReverse, SBV_ADDR_TYPE *pcData, SBV_ADDR_TYPE *cDataReverse, 
											u_int32_t nodeNum, SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE lDstAddr = lDataReverse[idx];
	if(lDstAddr == NULL_ADDR){
		return;
	}
	SBV *lDstSbv = &sbvPool[lDstAddr];
	if(lDstSbv->empty(sbvPool)){
		return;
	}


	SBV_ADDR_TYPE cUpdateAddr = cDataReverse[idx];
	u_int32_t lbuffer[MAX_EDGE_BUFFER_SIZE];
	int llen = lDstSbv->toArray(lbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < llen; ++i){

		auto lDst = lbuffer[i];
		SBV_ADDR_TYPE uDstAddr = uData[lDst];
		u_int32_t ubuffer[MAX_EDGE_BUFFER_SIZE];
		int ulen = sbvPool[uDstAddr].toArray(ubuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

		for(int j = 0; j < ulen; ++j){
			auto uDst = ubuffer[j];
			SBV_ADDR_TYPE pcDstAddr = pcData[uDst];
			if(sbvPool[pcDstAddr].test(idx, pcDstAddr, sbvPool)){
				if(sbvPool[cUpdateAddr].set(uDst, cUpdateAddr, sbvPool, sbvPoolIndex)){
					// printf("Adding reverse c-edge (%d,%d)\n", idx, uDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
	return;
}

__global__ void rule3Reverse(SBV_ADDR_TYPE *dDataReverse, SBV_ADDR_TYPE *sData, SBV_ADDR_TYPE *pcData, SBV_ADDR_TYPE *cDataReverse, 
					u_int32_t nodeNum, SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE dDstAddr = dDataReverse[idx];
	if(dDstAddr == NULL_ADDR){
		return;
	}

	SBV *dDstSbv = &sbvPool[dDstAddr];
	if(dDstSbv->empty(sbvPool)){
		return;
	}

	SBV_ADDR_TYPE cUpdateAddr = cDataReverse[idx];
	u_int32_t dbuffer[MAX_EDGE_BUFFER_SIZE];
	int dlen = dDstSbv->toArray(dbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < dlen; ++i){
		auto dDst = dbuffer[i];
		SBV_ADDR_TYPE sDstAddr = sData[dDst];
		u_int32_t sbuffer[MAX_EDGE_BUFFER_SIZE];
		int slen = sbvPool[sDstAddr].toArray(sbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

		for(int j = 0; j < slen; ++j){
			auto sDst = sbuffer[j];

			SBV_ADDR_TYPE pcDstAddr = pcData[sDst];
			if(sbvPool[pcDstAddr].test(idx, pcDstAddr, sbvPool)){
				if(sbvPool[cUpdateAddr].set(sDst, cUpdateAddr, sbvPool, sbvPoolIndex)){
					// printf("Adding reverse c-edge (%d,%d)\n", idx, sDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
}

__global__ void rule4Reverse(SBV_ADDR_TYPE *pdDataReverse, SBV_ADDR_TYPE *pData, SBV_ADDR_TYPE *dDataReverse, u_int32_t nodeNum, 
					SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed, const u_int32_t* nodeIdToMrId,
    		const u_int32_t* derivedIds, const u_int32_t* derivedOffsets, const u_int32_t* derivedCounts){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}
	
	SBV_ADDR_TYPE pdDstAddr = pdDataReverse[idx];
	if(pdDstAddr == NULL_ADDR){
		return;
	}
	SBV *pdDstSbv = &sbvPool[pdDstAddr];
	
	if(pdDstSbv->empty(sbvPool)){
		return;
	}

	SBV_ADDR_TYPE dDstAddr = dDataReverse[idx];
	SBV *dDstSbv = &sbvPool[dDstAddr];
	u_int32_t pdbuffer[MAX_EDGE_BUFFER_SIZE];
	int pdlen = pdDstSbv->toArray(pdbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	

	for(int i = 0; i < pdlen; ++i){
		auto pdDst = pdbuffer[i];
		SBV_ADDR_TYPE pDstAddr = pData[pdDst];
		u_int32_t pbuffer[MAX_EDGE_BUFFER_SIZE];
		int plen = sbvPool[pDstAddr].toArray(pbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
		for(int j = 0; j < plen; ++j){
			auto pDst = pbuffer[j];

			if(isaVersionofPtr(pDst, idx, nodeIdToMrId, nodeNum, derivedIds, derivedOffsets, derivedCounts)){
				if(dDstSbv->set(pdDst, dDstAddr, sbvPool, sbvPoolIndex)){
					// printf("Adding reverse d-edge (%d,%d)\n", idx, pdDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
	return;
}

__global__ void rule5Reverse(SBV_ADDR_TYPE *puData, SBV_ADDR_TYPE *pData, SBV_ADDR_TYPE *uData, u_int32_t nodeNum, SBV *sbvPool, 
					u_int32_t *sbvPoolIndex, int *changed, const u_int32_t* nodeIdToMrId,
    			const u_int32_t* derivedIds, const u_int32_t* derivedOffsets, const u_int32_t* derivedCounts){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}


	SBV_ADDR_TYPE puDstAddr = puData[idx];
	SBV_ADDR_TYPE pDstAddr = pData[idx];
	if(puDstAddr == NULL_ADDR || pDstAddr == NULL_ADDR){
		return;
	}

	SBV *puDstSbv = &sbvPool[puDstAddr];
	SBV *pDstSbv = &sbvPool[pDstAddr];

	if(puDstSbv->empty(sbvPool) || pDstSbv->empty(sbvPool)){
		return;
	}

	SBV_ADDR_TYPE uDstAddr = uData[idx];
	u_int32_t pubuffer[MAX_EDGE_BUFFER_SIZE];
	int pulen = puDstSbv->toArray(pubuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	u_int32_t pbuffer[MAX_EDGE_BUFFER_SIZE];
	int plen = pDstSbv->toArray(pbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < pulen; ++i){
		auto puDst = pubuffer[i];
		for(int j = 0; j < plen; ++j){
			auto pDst = pbuffer[j];
			// printf("(%d,%d,%d)\n", pDst, puDst, isaVersionofPtr(pDst, puDst, nodeIdToMrId, nodeNum, derivedIds, derivedOffsets, derivedCounts));
			if(isaVersionofPtr(pDst, puDst, nodeIdToMrId, nodeNum, derivedIds, derivedOffsets, derivedCounts)){
				if(sbvPool[uDstAddr].set(puDst, uDstAddr, sbvPool, sbvPoolIndex)){
					// printf("Adding u-edge (%d,%d)\n", idx, puDst);
					atomicExch(changed, 1);
				}
			}
		}
	}
}

__global__ void rule6Reverse(SBV_ADDR_TYPE *dDataReverse, SBV_ADDR_TYPE *kDataReverse, SBV_ADDR_TYPE *aData, u_int32_t nodeNum, 
					SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}
	
	SBV_ADDR_TYPE kDstAddr = kDataReverse[idx];
	SBV_ADDR_TYPE aDstAddr = aData[idx];
	// printf("kdstsbv %d", kDstSbv->empty());

	SBV *kDstSbv = &sbvPool[kDstAddr];
	if(kDstSbv->empty(sbvPool)){
		return;
	}

	u_int32_t kbuffer[MAX_EDGE_BUFFER_SIZE];
	int klen = kDstSbv->toArray(kbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	for(int i = 0; i < klen; ++i){
		auto kDst = kbuffer[i];

		SBV_ADDR_TYPE dDstAddr = dDataReverse[kDst];
		SBV *dDstSbv = &sbvPool[dDstAddr];
		// printf("dDstSbv %d", dDstSbv->empty());

		if(!dDstSbv->empty(sbvPool)){
			if(sbvPool[aDstAddr].set(kDst, aDstAddr, sbvPool, sbvPoolIndex)){
				// printf("Adding a-edge (%d,%d)\n", idx, kDst);
				atomicExch(changed, 1);
			}
		}
	}
}

__global__ void rule7Reverse(SBV_ADDR_TYPE *nDataReverse, SBV_ADDR_TYPE *kData, SBV_ADDR_TYPE *aData, SBV_ADDR_TYPE *cDataReverse, 
					u_int32_t nodeNum, SBV *sbvPool, u_int32_t *sbvPoolIndex, int *changed, const u_int32_t* nodeIdToMrId){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx == 0 || idx > nodeNum){
	  return;
	}

	SBV_ADDR_TYPE nDstAddr = nDataReverse[idx];
	SBV_ADDR_TYPE kDstAddr = kData[idx];

	SBV *nDstSbv = &sbvPool[nDstAddr];
	SBV *kDstSbv = &sbvPool[kDstAddr];
	if(nDstSbv->empty(sbvPool) || kDstSbv->empty(sbvPool)){
		return;
	}
	SBV_ADDR_TYPE cDstAddr = cDataReverse[idx];
	const u_int32_t selfMrId = getMrIdOfPtgNode(idx, nodeIdToMrId, nodeNum);

	u_int32_t nbuffer[MAX_EDGE_BUFFER_SIZE];
	int nlen = nDstSbv->toArray(nbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
	u_int32_t kbuffer[MAX_EDGE_BUFFER_SIZE];
	int klen = kDstSbv->toArray(kbuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);

	for(int i = 0; i < nlen; ++i){
		auto nDst = nbuffer[i];

		for(int j = 0; j < klen; ++j){
			auto kDst = kbuffer[j];

			SBV_ADDR_TYPE aDstAddr = aData[kDst];
			u_int32_t abuffer[MAX_EDGE_BUFFER_SIZE];
			int alen = sbvPool[aDstAddr].toArray(abuffer, sbvPool, MAX_EDGE_BUFFER_SIZE);
			for(int k = 0; k < alen; ++k){
				auto aDst = abuffer[k];
				if(getMrIdOfPtgNode(nDst, nodeIdToMrId, nodeNum) != getMrIdOfPtgNode(aDst, nodeIdToMrId, nodeNum)){
					if(sbvPool[cDstAddr].set(nDst, cDstAddr, sbvPool, sbvPoolIndex)){
						// printf("Adding reverse c-edge (%d,%d)\n", idx, nDst);
						atomicExch(changed, 1);
					}
				}
			}
		}
	}
}


void solveReverse(std::map<std::string, SBV_ADDR_TYPE*> dataMap, u_int32_t nodeNum, SBV *deviceSbvPool, u_int32_t *deviceSbvPoolIndex,
						const u_int32_t* nodeIdToMrId, const u_int32_t* derivedIds, const u_int32_t* derivedOffsets, 
						const u_int32_t* derivedCounts, bool debug = false, size_t threadsPerBlock = 16){

	DEBUG_LOG("Starting solving with reverse edges\n");

	int hostChanged = 1;
	int *deviceChanged;
	cudaMalloc(&deviceChanged, sizeof(int));

	size_t round = 0;
	// valid node id start from 1.
	size_t blockPerGrid = (nodeNum+threadsPerBlock) / threadsPerBlock;
	std::cout << blockPerGrid << " " << threadsPerBlock << "\n";
	while(hostChanged){
		++round;
		// DEBUG_LOG("Round: " << round);
		
		hostChanged = 0;
		CUDA_CHECK(cudaMemcpy(deviceChanged, &hostChanged, sizeof(int), cudaMemcpyHostToDevice));

		rule1Reverse<<<blockPerGrid,threadsPerBlock>>>(dataMap["p"], dataMap["c-1"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		cudaError_t err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
				std::cerr << "CUDA kernel error1: " << cudaGetErrorString(err) << std::endl;
		}
		rule2Reverse<<<blockPerGrid,threadsPerBlock>>>(dataMap["u"], dataMap["l-1"], dataMap["p_c"], dataMap["c-1"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
				std::cerr << "CUDA kernel error2: " << cudaGetErrorString(err) << std::endl;
		}
		rule3Reverse<<<blockPerGrid,threadsPerBlock>>>(dataMap["d-1"], dataMap["s"], dataMap["p_c"], dataMap["c-1"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
				std::cerr << "CUDA kernel error3: " << cudaGetErrorString(err) << std::endl;
		}
		rule4Reverse<<<blockPerGrid,threadsPerBlock>>>(dataMap["p_d-1"], dataMap["p"], dataMap["d-1"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged, nodeIdToMrId, derivedIds,
										derivedOffsets, derivedCounts);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
				std::cerr << "CUDA kernel error4: " << cudaGetErrorString(err) << std::endl;
		}
		rule5Reverse<<<blockPerGrid,threadsPerBlock>>>(dataMap["p_u"], dataMap["p"], dataMap["u"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged, nodeIdToMrId, derivedIds,
										derivedOffsets, derivedCounts);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
				std::cerr << "CUDA kernel error5: " << cudaGetErrorString(err) << std::endl;
		}
		rule6Reverse<<<blockPerGrid,threadsPerBlock>>>(dataMap["d-1"], dataMap["k-1"], dataMap["a"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
				std::cerr << "CUDA kernel error6: " << cudaGetErrorString(err) << std::endl;
		}
		rule7Reverse<<<blockPerGrid,threadsPerBlock>>>(dataMap["n-1"], dataMap["k"], dataMap["a"], dataMap["c-1"], nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceChanged, nodeIdToMrId);

		err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
				std::cerr << "CUDA kernel error7: " << cudaGetErrorString(err) << std::endl;
		}
    CUDA_CHECK(cudaMemcpy(&hostChanged, deviceChanged, sizeof(int), cudaMemcpyDeviceToHost));
	}
	
	DEBUG_LOG("Round: " << round);
	

}


void verifyResult(std::map<std::string, std::map<u_int32_t, std::set<u_int32_t>>> &edgeTypeToSbvGroupMap, std::string expectedOutputFileName){

	std::ifstream inFile(expectedOutputFileName+"-reversed.out");
	if(!inFile){
		std::cerr << "Cannot open file " << expectedOutputFileName << ".out\n";
		std::terminate();
	}

	std::set<std::tuple<u_int32_t, u_int32_t, std::string>> data;
	std::string line;
	while(std::getline(inFile, line)){
		std::stringstream ss(line);
		std::string element;
		std::vector<std::string> elements;
		while(std::getline(ss, element, ',')){
			elements.push_back(element);
		}
		u_int32_t from = static_cast<uint32_t>(std::stoul(elements[0]));
		u_int32_t to = static_cast<uint32_t>(std::stoul(elements[1]));
		std::tuple<u_int32_t, u_int32_t, std::string> t(from, to, elements[2]);
		data.insert(t);
	}

	inFile.close();

	std::set<std::tuple<u_int32_t, u_int32_t, std::string>> extra;
	// std::vector<std::string> edgeType{"p", "c", "s", "l", "p_d", "p_u", "p_c", "d", "u", "n", "k", "a"};
	// reverse edges
	std::set<std::string> edgeType{"p", "u", "s", "p_c", "p_u", "u", "a", "k", "c-1", "l-1", "d-1", "p_d-1", "n-1"};
	// size_t numberOfTotalEdges = 0;
	for(auto et : edgeType){
		if(!edgeTypeToSbvGroupMap.count(et)){
			continue;
		}
		auto &pts = edgeTypeToSbvGroupMap.at(et);
		for(const auto &[from, tos]: pts){
			for(auto to : tos){
				// check existance for edge from - et -> to in data
				std::tuple<u_int32_t, u_int32_t, std::string> t(from, to, et);
				// ++numberOfTotalEdges;
				auto iter = data.find(t);
				if(iter != data.end()){
          data.erase(iter);
        }
        else{
          extra.insert(t);
        }
			}
		}
	}


	bool error = false;
  if(!data.empty()){
    std::cerr << "ptgSbv missing following points-to entries:\n";
    for(auto [from, to, type] : data){
			// if(type != "p"){
			// 	continue;
			// }
      std::string fromName;
      std::stringstream ss(fromName); 
      ss << from;
      ss << " == " << type << " ==> ";
      ss << to;
      std::cerr << ss.str() << "\n";
    }
    error = true;
  }

  if(!extra.empty()){
    // std::cerr << "ptgSbv contains following extra entries:\n";
    // for(auto [from, to, type] : extra){
    //   std::string fromName;
    //   std::stringstream ss(fromName); 
    //   ss << from;
    //   ss << " == " << type << " ==> ";
    //   ss << to;
    //   std::cerr << ss.str() << "\n";
    // }
    error = true;
  }

  if(error){
    std::cerr << "Fail to verify. " << data.size() << " edges are missing, " <<  extra.size() << " extra edges\n";
  }
  else{
    std::cout << "Verify succeed.\n";
  }

	// std::cout << "Total contains " << numberOfTotalEdges << " edges\n";

}


int gpamain(const std::string ptgFileName){

	auto ptg = readPtg(ptgFileName);
	std::cout << "Initial ptg contains " << ptg.size << " edges\n";
	// auto edgeTypeToSbvGroupMap = preprocessPtgForGpu(ptg);	
	auto edgeTypeToSbvGroupMap = preprocessPtgForGpuWithReverseEdges(ptg);	

	auto start = std::chrono::high_resolution_clock::now();
	std::cout << "1\n";

	// set up gpu memory pool.
	SBV* deviceSbvPool;
	CUDA_CHECK(cudaMalloc(&deviceSbvPool, sizeof(SBV) * MAX_ALLOCATED_SBV_NUM));
	u_int32_t index = 0;
	std::cout << "2\n";
	
	// copy initialized SBVs from CPU to GPU.
	std::map<std::string, SBV_ADDR_TYPE*> dataMap;
	std::set<std::string> edgeType{"p", "u", "s", "p_c", "p_u", "u", "a", "k", "c-1", "l-1", "d-1", "p_d-1", "k-1", "n-1"};
	// node 0 is reserved for nullptr. Valid node index starts from 1.
	size_t allocatedCount = ptg.nodeNum+1;
	std::cout << "3\n";

	size_t used = 0;
	for(auto et : edgeType){
		const auto& hostHeads = edgeTypeToSbvGroupMap[et];

		size_t totalBlocks = 0;
    for(size_t i = 0; i < allocatedCount; ++i){
			for(auto* p = hostHeads[i]; p; p = p->next){
				++totalBlocks;
			}
		}
         
		std::vector<SBV> h_blocks; 
		h_blocks.reserve(totalBlocks);
    std::vector<SBV_ADDR_TYPE> head_idx(allocatedCount, -1);

		for(size_t i = 0; i < allocatedCount; ++i) {
			auto* p = hostHeads[i];
			SBV_ADDR_TYPE prev = -1;
			while(p){
				SBV node;
				node.base = p->base;
				for(int w = 0; w < 29; ++w){
					node.bits[w] = p->bits[w];
				} 
				node.next = -1;

				SBV_ADDR_TYPE idx = static_cast<SBV_ADDR_TYPE>(h_blocks.size());
				h_blocks.emplace_back(node);
				if(prev == -1){
					head_idx[i] = idx;
				} 
				else{
					h_blocks[prev].next = idx;
				}           
				prev = idx; 
				p = p->next;
			}
    }

		for (auto& n : h_blocks) {
			if(n.next != -1){
				n.next += static_cast<SBV_ADDR_TYPE>(used);
			} 
    }
    for (auto& h : head_idx) {
			if(h != -1){
				h += static_cast<SBV_ADDR_TYPE>(used);
			} 
    }

		if(!h_blocks.empty()) {
      CUDA_CHECK(cudaMemcpy(deviceSbvPool + used, h_blocks.data(), h_blocks.size() * sizeof(SBV), cudaMemcpyHostToDevice));
    }

    // std::vector<SBV*> h_head_ptrs(allocatedCount, nullptr);
    // for (size_t i = 0; i < allocatedCount; ++i) {
    //     if (h_heads[i] != -1) h_head_ptrs[i] = deviceSbvPool + h_heads[i];
    // }

    SBV_ADDR_TYPE* d_head_idx = nullptr;  
    CUDA_CHECK(cudaMalloc(&d_head_idx, allocatedCount * sizeof(SBV_ADDR_TYPE)));
    CUDA_CHECK(cudaMemcpy(d_head_idx, head_idx.data(), allocatedCount * sizeof(SBV_ADDR_TYPE), cudaMemcpyHostToDevice));

    dataMap[et] = d_head_idx;

    // Advance global offset
    used += h_blocks.size();

	}

	// CUDA_CHECK(cudaMemcpy(d_sbvPoolIndex, &used, sizeof(uint32_t), cudaMemcpyHostToDevice));







	std::cout << "4\n";

	u_int32_t* deviceSbvPoolIndex;
	CUDA_CHECK(cudaMalloc(&deviceSbvPoolIndex, sizeof(u_int32_t)));
	// index refers to the address of the first available SBV block in the GPU pool after initialization.
	CUDA_CHECK(cudaMemcpy(deviceSbvPoolIndex, &used, sizeof(uint32_t), cudaMemcpyHostToDevice));

	std::vector<u_int32_t> hostNodeIdToMrId(allocatedCount, 0);
	for (const auto &[nodeId, tup] : ptg.nodeId2mrIdAndVerMap){
    hostNodeIdToMrId[nodeId] = std::get<1>(tup);  // the mrId
	}

	u_int32_t* deviceNodeIdToMrId;
	CUDA_CHECK(cudaMalloc(&deviceNodeIdToMrId, hostNodeIdToMrId.size() * sizeof(u_int32_t)));
	CUDA_CHECK(cudaMemcpy(deviceNodeIdToMrId, hostNodeIdToMrId.data(), allocatedCount * sizeof(u_int32_t), cudaMemcpyHostToDevice));

	std::cout << "5\n";

	std::vector<u_int32_t> hostDerivedIds, hostDerivedOffsets, hostDerivedCounts;

	for (u_int32_t i = 0; i < allocatedCount; ++i){
			auto it = ptg.ptrIdToDerivedMrIdsMap.find(i);
			if(it == ptg.ptrIdToDerivedMrIdsMap.end()){
					hostDerivedOffsets.push_back(hostDerivedIds.size());
					hostDerivedCounts.push_back(0);
			} 
			else{
				hostDerivedOffsets.push_back(hostDerivedIds.size());
				hostDerivedCounts.push_back(it->second.size());
				for(u_int32_t id : it->second) {
						hostDerivedIds.push_back(id);
				}
			}
	}

	std::cout << "6\n";
	u_int32_t *deviceDerivedIds, *deviceDerivedOffsets, *deviceDerivedCounts;

	CUDA_CHECK(cudaMalloc(&deviceDerivedIds, hostDerivedIds.size() * sizeof(u_int32_t)));
	CUDA_CHECK(cudaMalloc(&deviceDerivedOffsets, hostDerivedOffsets.size() * sizeof(u_int32_t)));
	CUDA_CHECK(cudaMalloc(&deviceDerivedCounts, hostDerivedCounts.size() * sizeof(u_int32_t)));
	CUDA_CHECK(cudaMemcpy(deviceDerivedIds, hostDerivedIds.data(), hostDerivedIds.size() * sizeof(u_int32_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(deviceDerivedOffsets, hostDerivedOffsets.data(), hostDerivedOffsets.size() * sizeof(u_int32_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(deviceDerivedCounts, hostDerivedCounts.data(), hostDerivedCounts.size() * sizeof(u_int32_t), cudaMemcpyHostToDevice));

	cudaDeviceSynchronize();



	auto preprocessed = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(preprocessed - start);
  std::cout << "Preprocessing time: " << duration.count() << " ms\n";

	cudaEvent_t cudastart, cudastop;
	cudaEventCreate(&cudastart);
	cudaEventCreate(&cudastop);

	cudaEventRecord(cudastart);

	// solve(dataMap, ptg.nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceNodeIdToMrId, deviceDerivedIds, deviceDerivedOffsets,
	// 			deviceDerivedCounts, true, 1);

	solveReverse(dataMap, ptg.nodeNum, deviceSbvPool, deviceSbvPoolIndex, deviceNodeIdToMrId, deviceDerivedIds, deviceDerivedOffsets,
				deviceDerivedCounts, true);

	cudaEventRecord(cudastop);

	cudaEventSynchronize(cudastop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, cudastart, cudastop);

	std::cout << "Solving time (GPU only): " << milliseconds << " ms\n";

	cudaEventDestroy(cudastart);
	cudaEventDestroy(cudastop);

	cudaDeviceSynchronize();
	auto solved = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(solved - preprocessed);
  std::cout << "Solving time: " << duration.count() << " ms\n";

	auto end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Total time: " << duration.count() << " ms\n";

	u_int32_t totalSbvBlockNum;
	CUDA_CHECK(cudaMemcpy(&totalSbvBlockNum, deviceSbvPoolIndex, sizeof(u_int32_t), cudaMemcpyDeviceToHost));

	std::cout << "Total number of sbv block allocated is " << totalSbvBlockNum << "\n";








	std::map<std::string, std::map<u_int32_t, std::set<u_int32_t>>> edgeResultMap;
	std::vector<SBV> hostPool;
	hostPool.resize(MAX_ALLOCATED_SBV_NUM);
	CUDA_CHECK(cudaMemcpy(hostPool.data(), deviceSbvPool, MAX_ALLOCATED_SBV_NUM * sizeof(SBV), cudaMemcpyDeviceToHost));



	// for(const std::string& label : edgeType) {
  //   SBV_ADDR_TYPE *deviceLabelData = dataMap[label];

  //   // Host array of device pointers
  //   std::vector<SBV_ADDR_TYPE> hostSbvPtrArray(allocatedCount);
  //   CUDA_CHECK(cudaMemcpy(hostSbvPtrArray.data(), deviceLabelData, allocatedCount * sizeof(SBV_ADDR_TYPE), cudaMemcpyDeviceToHost));

  //   for(u_int32_t src = 0; src < allocatedCount; ++src){
	// 		SBV_ADDR_TYPE deviceAddr = hostSbvPtrArray[src];
	// 		if(deviceAddr == NULL_ADDR){
	// 			continue;
	// 		} 

			// Copy chain of SBVs to host
			// std::vector<SBV> chain;
			// SBV curr;
			// CUDA_CHECK(cudaMemcpy(&curr, deviceSbv, sizeof(SBV), cudaMemcpyDeviceToHost));
			// chain.push_back(curr);
			// while (curr.next) {
			// 		CUDA_CHECK(cudaMemcpy(&curr, curr.next, sizeof(SBV), cudaMemcpyDeviceToHost));
			// 		chain.push_back(curr);
			// }

			// Decode all IDs
	// 		SBV_ADDR_TYPE curr = deviceAddr;
	// 		while(curr != NULL_ADDR){
	// 			SBV &sbv = hostPool[curr];
	// 			for(int i = 0; i < WORDS_PER_BLOCK; ++i){
	// 				u_int32_t word = sbv.bits[i];
	// 				if(!word){
	// 					continue;
	// 				} 
	// 				for(int bit = 0; bit < 32; ++bit){
	// 					if(word & (1u << bit)){
	// 						u_int32_t baseId = sbv.base;
	// 						u_int32_t wordId = 28 - i;
	// 						u_int32_t inWordId = bit;
	// 						u_int32_t id = baseId * 928 + wordId * 32 + inWordId;
	// 						edgeResultMap[label][src].insert(id);
	// 					}
	// 				}
	// 			}

	// 			curr = sbv.next;
	// 		}
  //   }
	// }

	// verifyResult(edgeResultMap, ptgFileName);

	// for(const auto &[label, pts] : edgeResultMap){
	// 	if(label != "p"){
	// 		continue;
	// 	}
  //   for(const auto &[src, dsts] : pts){
	// 		if(dsts.empty()){
	// 			continue;
	// 		}
	// 		std::cout << src << " ==> {";
	// 		for(auto dst : dsts){
	// 			std::cout << dst << " ";
	// 		}
	// 		std::cout << "}\n";
  //   }
	// }

	// clean up

	std::cout << "7\n";
	CUDA_CHECK(cudaFree(deviceSbvPool));
	CUDA_CHECK(cudaFree(deviceSbvPoolIndex));
	CUDA_CHECK(cudaFree(deviceNodeIdToMrId));
	CUDA_CHECK(cudaFree(deviceDerivedIds));
	CUDA_CHECK(cudaFree(deviceDerivedOffsets));
	CUDA_CHECK(cudaFree(deviceDerivedCounts));

	std::cout << "8\n";


	// for(const auto &[et, ptrs] : dataMap){
	// 	CUDA_CHECK(cudaFree(ptrs));
	// }
	
	std::cout << "9\n";


	for(auto et : edgeType){
		for(auto sbv : edgeTypeToSbvGroupMap[et]){
			delete sbv;
		}
	}

	std::cout << "10\n";


	return 0;
	
}