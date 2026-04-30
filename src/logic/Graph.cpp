#include "Graph.hpp"
#include <limits>
#include <tuple>
#include <iostream>
#include <geometry/FruchtermanReingold_custom.h>
#include <vector>

void Graph::addNode(int id, double x, double y) {
    nodes.push_back({id, x, y}); // push_back adds an item to the end of a vector
}

void Graph::addEdge(int id, int u_index, int v_index) {
    edges.push_back({id, u_index, v_index, -1});
}

Graph Graph::createPlanarizedGraph() const {
    Graph newGraph;
    // ... logic for planarization ...
    return newGraph;
}

std::tuple<double, double, double, double> Graph::getBounds() const {
    if (nodes.empty()) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& node : nodes) {
        if (node.x < minX) minX = node.x;
        if (node.x > maxX) maxX = node.x;
        if (node.y < minY) minY = node.y;
        if (node.y > maxY) maxY = node.y;
    }

    return {minX, minY, maxX, maxY};
}



void Graph::applyInitialForceDirectedLayout(int iterations) {
    if (nodes.empty()) return;

    std::cout << "[INIT] Starting Fruchterman-Reingold layout for " << iterations << " iterations." << std::endl;

    std::tuple<double, double, double, double> bounds = getBounds();

    double minX = std::get<0>(bounds);
    double minY = std::get<1>(bounds);
    double maxX = std::get<2>(bounds);
    double maxY = std::get<3>(bounds);

    FruchtermanReingold fr;
    fr.initialize(maxX - minX, minX, minY);

    std::vector<FruchtermanReingold::Vector2D> repForces;
    std::vector<FruchtermanReingold::Vector2D> attrForces;

    double initialTemperature = 100.0; 

    for (int i = 0; i < iterations; ++i) {
        // Temperature cools down linearly as iterations progress
        double temp = initialTemperature * (1.0 - static_cast<double>(i) / iterations);
        if (temp < 0.1) temp = 0.1;

        // Note: If m_pGraph (PlanarizedGraph) is not strictly compatible with 'Graph', 
        // you may need to map m_pGraph's nodes to your FruchtermanReingold graph type here.
        
        // 1. Calculate Forces
        // If your graph is huge, use calculateApproximateRepulsiveForces with your grid
        fr.calculateExactRepulsiveForces(*this, repForces); 
        fr.calculateAttractiveForces(*this, attrForces);

        // 2. Apply Displacement
        fr.updateNodePositions(*this, repForces, attrForces, temp);

        std::cout << "[INIT] Force-directed layout complete." << std::endl;
    }
}