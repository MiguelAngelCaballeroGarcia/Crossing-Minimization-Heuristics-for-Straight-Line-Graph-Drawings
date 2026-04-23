#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>
#include "Graph.hpp" // Your original graph
#include "SpatialGrid.hpp" // For spatial indexing and incremental updates

class PlanarizedGraph {
private:
    SpatialGrid grid;

    int nextNodeId = 0;
    int nextEdgeId = 0;
    std::vector<int> freeNodeIds;
    std::vector<int> freeEdgeIds;

    // Helper to safely swap a neighbor pointer in either an ORIGINAL or CROSSING node
    void updateNodeNeighbor(int nodeId, int oldNeighbor, int newNeighbor, int origEdgeId);

    // Helpers to manage the PlanarEdges safely
    void removePlanarEdge(int u_id, int v_id, int origEdgeId);
    void addPlanarEdge(int u_id, int v_id, int origEdgeId);

    // Math helper for createCrossing
    bool isPointOnSegment(double px, double py, int u_id, int v_id);

    void ensureNodeCapacity(int id);
    void ensureEdgeCapacity(int id);
    void ensureOriginalEdgeCapacity(int id);
    void deactivateNode(int id);
    void deactivateEdge(int id);
public:
    struct IntersectionData {
        double x, y;
        int edge1_id;
        int edge2_id;
        double t1; // Position along edge1 (0.0 at start, 1.0 at end)
        double t2; // Position along edge2 (0.0 at start, 1.0 at end)
    };

    const SpatialGrid& getGrid() const;

    PlanarizedGraph(const Graph& originalGraph, const std::vector<IntersectionData>& intersections);

    int getNextEdgeId();
    int getNextNodeId();

    enum class NodeType { ORIGINAL, CROSSING };

    struct PlanarNode {
        int id;
        double x, y;
        NodeType type;
        std::vector<int> incident_planar_edges;
        
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
    std::vector<std::vector<int>> originalEdgeToPlanarEdges;
    std::unordered_map<std::uint64_t, int> crossingPairToNode;

    std::vector<PlanarNode> nodes;
    std::vector<std::uint8_t> nodeActive;
    std::vector<PlanarEdge> edges;
    std::vector<std::uint8_t> edgeActive;

    bool hasNode(int nodeId) const;
    bool hasEdge(int edgeId) const;
    PlanarNode& getNode(int nodeId);
    const PlanarNode& getNode(int nodeId) const;
    PlanarEdge& getEdge(int edgeId);
    const PlanarEdge& getEdge(int edgeId) const;

    template <typename Func>
    void forEachNode(Func&& fn) {
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            if (nodeActive[i] == 0) continue;
            fn(i, nodes[i]);
        }
    }

    template <typename Func>
    void forEachNode(Func&& fn) const {
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            if (nodeActive[i] == 0) continue;
            fn(i, nodes[i]);
        }
    }

    template <typename Func>
    void forEachEdge(Func&& fn) {
        for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
            if (edgeActive[i] == 0) continue;
            fn(i, edges[i]);
        }
    }

    template <typename Func>
    void forEachEdge(Func&& fn) const {
        for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
            if (edgeActive[i] == 0) continue;
            fn(i, edges[i]);
        }
    }

    // Helper to find where a crossing sits between two nodes on an edge
    void insertCrossingIntoEdgeChain(int edgeId, int crossingNodeIdx);
    
    void destroyCrossing(int crossingNodeIdx);
    void createCrossing(int edge1Id, int edge2Id, double x, double y);
    int getCrossingNodeForPair(int edgeA, int edgeB) const;

    // Moves a node and updates grid occupancy for both node and incident planar edges.
    void updateNodePosition(int nodeId, double newX, double newY);

    // Count total crossing nodes in the graph
    int countTotalCrossings() const;
};