#include "Graph.hpp"

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