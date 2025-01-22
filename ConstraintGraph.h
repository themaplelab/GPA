#ifndef GPU_FSPA_CONSTRAINT_GRAPH_H
#define GPU_FSPA_CONSTRAINT_GRAPH_H

#include "Constraints.h"

#include <vector>
#include <memory>
#include <set>
#include <map>




class ConstraintGraph{

    std::set<size_t> nodes;
    std::map<size_t, std::set<size_t>> pEdges;
    std::map<size_t, std::set<size_t>> cEdges;
    std::map<size_t, std::set<size_t>> sEdges;
    std::map<size_t, std::set<size_t>> lEdges;

    public:
        ConstraintGraph(std::vector<Constraint> &constraints){
            for(auto c : constraints){
                auto lhs = c.getLhs();
                auto rhs = c.getRhs();
                auto type = c.getType();

                if(!nodes.count(lhs)){
                    nodes.insert(lhs);
                }

                if(!nodes.count(rhs)){
                    nodes.insert(rhs);
                }

                switch (type){
                case Constraint::ConstraintType::PointsTo:
                    pEdges[lhs].insert(rhs);
                    break;

                case Constraint::ConstraintType::Copy:
                    cEdges[lhs].insert(rhs);
                    break;

                case Constraint::ConstraintType::Store:
                    sEdges[lhs].insert(rhs);
                    break;

                case Constraint::ConstraintType::Load:
                    lEdges[lhs].insert(rhs);
                    break;
                
                default:
                    break;
                }

            }
        }

        size_t getNodeNumbers() const {
            return nodes.size();
        }

        size_t getPedgeNumbers() const {
            size_t res = 0;
            for(auto p : pEdges){
                res += p.second.size();
            }

            return res;
        }

        size_t getCedgeNumbers() const {
            size_t res = 0;
            for(auto p : cEdges){
                res += p.second.size();
            }

            return res;
        }

        size_t getSedgeNumbers() const {
            size_t res = 0;
            for(auto p : sEdges){
                res += p.second.size();
            }

            return res;
        }

        size_t getLedgeNumbers() const {
            size_t res = 0;
            for(auto p : lEdges){
                res += p.second.size();
            }

            return res;
        }

        std::set<size_t> getNodes() const {
            return nodes;
        }

        std::map<size_t, std::set<size_t>> getPedges() const{
            return pEdges;
        }

        std::map<size_t, std::set<size_t>> getSedges() const{
            return sEdges;
        }

        std::map<size_t, std::set<size_t>> getLedges() const{
            return lEdges;
        }

        std::map<size_t, std::set<size_t>> getCedges() const{
            return cEdges;
        }

        bool addCedge(size_t from, size_t to){
            auto changed = cEdges[from].insert(to).second;
            return changed;
        }

        bool addPedge(size_t from, size_t to){
            auto changed = pEdges[from].insert(to).second;
            return changed;
        }

};


#endif