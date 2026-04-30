#include "FruchtermanReingold_custom.h"
#include "../logic/Graph.hpp"
#include <algorithm>
#include <cmath>

FruchtermanReingold::FruchtermanReingold()
    : boxLength(100.0),
      minX(0.0),
      minY(0.0),
      optimalDistance(1.0),
      repulsiveStrength(1000.0),
      attractiveStrength(0.1),
      maxDisplacement(10.0) {
}

void FruchtermanReingold::initialize(double bl, double mx, double my) {
    boxLength = bl;
    minX = mx;
    minY = my;
    // Optimization parameters (gridQuotient, etc.) are removed for simplicity
}

FruchtermanReingold::Vector2D FruchtermanReingold::calculateRepulsiveForce(const FruchtermanReingold::Vector2D& pos1, 
                                                                           const FruchtermanReingold::Vector2D& pos2) const {
    FruchtermanReingold::Vector2D delta = pos2 - pos1;
    double distance = delta.magnitude();
    
    if (distance < 1e-6) {
        return FruchtermanReingold::Vector2D(0, 0);
    }
    
    // F_rep = (k^2 / distance) * strength
    double repulsiveForce = (optimalDistance * optimalDistance * repulsiveStrength) / distance;
    
    return delta.normalize() * repulsiveForce;
}

FruchtermanReingold::Vector2D FruchtermanReingold::calculateAttractiveForce(const FruchtermanReingold::Vector2D& pos1, 
                                                                            const FruchtermanReingold::Vector2D& pos2) const {
    FruchtermanReingold::Vector2D delta = pos2 - pos1;
    double distance = delta.magnitude();
    
    if (distance < 1e-6) {
        return FruchtermanReingold::Vector2D(0, 0);
    }
    
    // F_attr = (distance^2 / k) * strength
    double attractiveForce = (distance * distance * attractiveStrength) / optimalDistance;
    
    // Pulls pos2 toward pos1 (negative delta)
    return delta.normalize() * (-attractiveForce);
}

void FruchtermanReingold::calculateExactRepulsiveForces(const Graph& graph, 
                                                        std::vector<FruchtermanReingold::Vector2D>& repulsiveForces) {
    int nodeCount = static_cast<int>(graph.nodes.size());
    repulsiveForces.assign(nodeCount, FruchtermanReingold::Vector2D(0, 0));
    
    // Exact O(N^2) pairwise repulsion
    for (int i = 0; i < nodeCount; ++i) {
        for (int j = i + 1; j < nodeCount; ++j) {
            FruchtermanReingold::Vector2D pos_i(graph.nodes[i].x, graph.nodes[i].y);
            FruchtermanReingold::Vector2D pos_j(graph.nodes[j].x, graph.nodes[j].y);
            
            FruchtermanReingold::Vector2D repulsive = calculateRepulsiveForce(pos_i, pos_j);
            
            repulsiveForces[j] = repulsiveForces[j] + repulsive;
            repulsiveForces[i] = repulsiveForces[i] + (repulsive * -1.0);
        }
    }
}

void FruchtermanReingold::calculateAttractiveForces(const Graph& graph,
                                                    std::vector<FruchtermanReingold::Vector2D>& attractiveForces) {
    int nodeCount = static_cast<int>(graph.nodes.size());
    attractiveForces.assign(nodeCount, FruchtermanReingold::Vector2D(0, 0));
    
    for (const auto& edge : graph.edges) {
        int u_idx = edge.u_id;
        int v_idx = edge.v_id;
        
        if (u_idx < 0 || u_idx >= nodeCount || v_idx < 0 || v_idx >= nodeCount) continue;
        
        FruchtermanReingold::Vector2D pos_u(graph.nodes[u_idx].x, graph.nodes[u_idx].y);
        FruchtermanReingold::Vector2D pos_v(graph.nodes[v_idx].x, graph.nodes[v_idx].y);
        
        FruchtermanReingold::Vector2D attractive = calculateAttractiveForce(pos_u, pos_v);
        
        attractiveForces[v_idx] = attractiveForces[v_idx] + attractive;
        attractiveForces[u_idx] = attractiveForces[u_idx] + (attractive * -1.0);
    }
}

void FruchtermanReingold::updateNodePositions(Graph& graph,
                                              const std::vector<FruchtermanReingold::Vector2D>& repulsiveForces,
                                              const std::vector<FruchtermanReingold::Vector2D>& attractiveForces,
                                              double temperature) {
    int nodeCount = static_cast<int>(graph.nodes.size());
    
    for (int i = 0; i < nodeCount; ++i) {
        FruchtermanReingold::Vector2D totalForce = repulsiveForces[i] + attractiveForces[i];
        double forceMagnitude = totalForce.magnitude();
        
        // Displacement is limited by the current "temperature"
        double displacement = std::min(forceMagnitude * temperature, maxDisplacement);
        
        if (forceMagnitude > 1e-6) {
            FruchtermanReingold::Vector2D direction = totalForce.normalize();
            graph.nodes[i].x += direction.x * displacement;
            graph.nodes[i].y += direction.y * displacement;
        }
    }
}