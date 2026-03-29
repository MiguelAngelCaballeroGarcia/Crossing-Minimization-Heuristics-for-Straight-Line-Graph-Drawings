#include "Graph.hpp"
#include <limits>
#include <tuple>

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