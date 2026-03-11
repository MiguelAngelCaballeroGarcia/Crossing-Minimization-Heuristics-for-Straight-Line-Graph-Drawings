#pragma once

#include <vector>
#include <unordered_map>
#include "Graph.hpp" // Your original graph
#include "SpatialGrid.hpp" // For spatial indexing and incremental updates

class PlanarizedGraph {
private:
    SpatialGrid grid;

    int nextNodeId = 0;
    int nextEdgeId = 0;

    // Helper to safely swap a neighbor pointer in either an ORIGINAL or CROSSING node
    void updateNodeNeighbor(int nodeId, int oldNeighbor, int newNeighbor, int origEdgeId);

    // Helpers to manage the PlanarEdges safely
    void removePlanarEdge(int u_id, int v_id, int origEdgeId);
    void addPlanarEdge(int u_id, int v_id, int origEdgeId);

    // Math helper for createCrossing
    bool isPointOnSegment(double px, double py, int u_id, int v_id);
public:
    struct IntersectionData {
        double x, y;
        int edge1_id;
        int edge2_id;
    };

    PlanarizedGraph(const Graph& originalGraph, const std::vector<IntersectionData>& intersections);

    int getNextEdgeId() { return nextEdgeId++; }
    int getNextNodeId() { return nextNodeId++; }

    enum class NodeType { ORIGINAL, CROSSING };

    struct PlanarNode {
        int id;
        double x, y;
        NodeType type;
        
        // For ORIGINAL nodes:
        int original_node_id = -1;
        std::vector<int> adjacent_planar_nodes; // Can be any length (Degree D)

        // For CROSSING nodes:
        // We explicitly pair the neighbors by their original edge
        int original_edge_1 = -1;
        int e1_neighbor_prev = -1; // The neighbor towards the 'u' of original edge 1
        int e1_neighbor_next = -1; // The neighbor towards the 'v' of original edge 1

        int original_edge_2 = -1;
        int e2_neighbor_prev = -1; 
        int e2_neighbor_next = -1;
    };

    struct PlanarEdge {
        int id;
        int u_id; // index in PlanarizedGraph::nodes
        int v_id; // index in PlanarizedGraph::nodes
        
        // EVERY planar edge must know which original edge it is a part of
        int original_edge_id; 
    };

    // Fast lookup: Original Edge ID -> List of Planar Edge IDs (its sub-segments)
    std::unordered_map<int, std::vector<int>> originalEdgeToPlanarEdges;

    std::unordered_map<int, PlanarNode> nodes;
    std::unordered_map<int, PlanarEdge> edges;

    // Helper to find where a crossing sits between two nodes on an edge
    void insertCrossingIntoEdgeChain(int edgeId, int crossingNodeIdx);
    
    // The two methods you requested
    void destroyCrossing(int crossingNodeIdx);
    void createCrossing(int edge1Id, int edge2Id, double x, double y);
};