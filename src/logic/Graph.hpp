#pragma once
#include <vector>

class Graph {
public:
    struct Node {
        int id;
        double x, y;
    };

    struct Edge {
        int id;
        int u_id;
        int v_id;
        int original_edge_id = -1; 
    };

    std::vector<Node> nodes;
    std::vector<Edge> edges;

    // Adds a node and returns its index
    void addNode(int id, double x, double y);

    // Adds an edge between two node indices
    void addEdge(int id, int u_index, int v_index);

    // A method that computes crossings, splits edges, and returns a totally new graph
    Graph createPlanarizedGraph() const; 
};