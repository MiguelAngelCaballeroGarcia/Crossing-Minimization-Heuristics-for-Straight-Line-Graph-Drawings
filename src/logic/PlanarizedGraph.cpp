#include "PlanarizedGraph.hpp"
#include <algorithm>

void PlanarizedGraph::destroyCrossing(int crossingNodeIdx) {
    PlanarNode& cNode = nodes[crossingNodeIdx];

    // --- HEAL EDGE 1 ---
    int e1 = cNode.original_edge_1;
    int prev1 = cNode.e1_neighbor_prev;
    int next1 = cNode.e1_neighbor_next;

    // 1. Tell the neighbors to point to each other instead of cNode
    updateNodeNeighbor(prev1, crossingNodeIdx, next1, e1);
    updateNodeNeighbor(next1, crossingNodeIdx, prev1, e1);

    // 2. Delete old sub-segments and create the new merged segment
    removePlanarEdge(prev1, crossingNodeIdx, e1);
    removePlanarEdge(crossingNodeIdx, next1, e1);
    addPlanarEdge(prev1, next1, e1);

    // --- HEAL EDGE 2 ---
    int e2 = cNode.original_edge_2;
    int prev2 = cNode.e2_neighbor_prev;
    int next2 = cNode.e2_neighbor_next;

    updateNodeNeighbor(prev2, crossingNodeIdx, next2, e2);
    updateNodeNeighbor(next2, crossingNodeIdx, prev2, e2);

    removePlanarEdge(prev2, crossingNodeIdx, e2);
    removePlanarEdge(crossingNodeIdx, next2, e2);
    addPlanarEdge(prev2, next2, e2);

    // --- CLEANUP ---
    nodes.erase(crossingNodeIdx); // O(1) removal
}

void PlanarizedGraph::createCrossing(int edge1Id, int edge2Id, double x, double y) {
    int cId = getNextNodeId();
    
    // Initialize the new CROSSING node
    PlanarNode cNode;
    cNode.id = cId; cNode.x = x; cNode.y = y;
    cNode.type = NodeType::CROSSING;
    cNode.original_edge_1 = edge1Id;
    cNode.original_edge_2 = edge2Id;

    int affectedEdges[2] = {edge1Id, edge2Id};

    for (int eId : affectedEdges) {
        // Find which sub-segment of 'eId' contains the point (x, y)
        int targetEdgeId = -1;
        int u = -1, v = -1;

        for (int pEdgeId : originalEdgeToPlanarEdges[eId]) {
            PlanarEdge& pe = edges[pEdgeId];
            if (isPointOnSegment(x, y, pe.u_id, pe.v_id)) {
                targetEdgeId = pEdgeId;
                u = pe.u_id;
                v = pe.v_id;
                break; 
            }
        }

        if (targetEdgeId == -1) continue; // Safety check (should mathematically never happen)

        // 1. Update the new crossing's pointers
        if (eId == edge1Id) {
            cNode.e1_neighbor_prev = u;
            cNode.e1_neighbor_next = v;
        } else {
            cNode.e2_neighbor_prev = u;
            cNode.e2_neighbor_next = v;
        }

        // 2. Tell the old neighbors to point to the new crossing
        updateNodeNeighbor(u, v, cId, eId);
        updateNodeNeighbor(v, u, cId, eId);

        // 3. Update the PlanarEdges
        removePlanarEdge(u, v, eId);
        addPlanarEdge(u, cId, eId);
        addPlanarEdge(cId, v, eId);
    }

    // Save the new node to the map
    nodes[cId] = cNode;
}

void PlanarizedGraph::updateNodeNeighbor(int nodeId, int oldNeighbor, int newNeighbor, int origEdgeId) {
    PlanarNode& node = nodes[nodeId];

    if (node.type == NodeType::ORIGINAL) {
        // For original nodes, we just replace the ID in the adjacency list
        for (int& adj : node.adjacent_planar_nodes) {
            if (adj == oldNeighbor) {
                adj = newNeighbor;
                break;
            }
        }
    } else {
        // For crossing nodes, we must be careful to update the correct original edge track
        if (node.original_edge_1 == origEdgeId) {
            if (node.e1_neighbor_prev == oldNeighbor) node.e1_neighbor_prev = newNeighbor;
            else if (node.e1_neighbor_next == oldNeighbor) node.e1_neighbor_next = newNeighbor;
        }
        if (node.original_edge_2 == origEdgeId) {
            if (node.e2_neighbor_prev == oldNeighbor) node.e2_neighbor_prev = newNeighbor;
            else if (node.e2_neighbor_next == oldNeighbor) node.e2_neighbor_next = newNeighbor;
        }
    }
}

void PlanarizedGraph::removePlanarEdge(int u_id, int v_id, int origEdgeId) {
    auto& subEdges = originalEdgeToPlanarEdges[origEdgeId];
    for (auto it = subEdges.begin(); it != subEdges.end(); ++it) {
        PlanarEdge& pe = edges[*it];
        // Order doesn't matter for an undirected edge
        if ((pe.u_id == u_id && pe.v_id == v_id) || (pe.u_id == v_id && pe.v_id == u_id)) {
            edges.erase(*it);
            subEdges.erase(it);
            break; // Stop after deleting to not invalidate the iterator
        }
    }
}

void PlanarizedGraph::addPlanarEdge(int u_id, int v_id, int origEdgeId) {
    int newId = getNextEdgeId();
    edges[newId] = {newId, u_id, v_id, origEdgeId};
    originalEdgeToPlanarEdges[origEdgeId].push_back(newId);
}

// Math helper: Checks if (px, py) lies within the bounding box of segment u-v
bool PlanarizedGraph::isPointOnSegment(double px, double py, int u_id, int v_id) {
    PlanarNode& u = nodes[u_id];
    PlanarNode& v = nodes[v_id];
    
    // We add a tiny epsilon (1e-7) because floating point math is never perfect.
    // Since (px, py) is an intersection point, it is mathematically guaranteed to be on the line,
    // so we only need to check if it's "between" the endpoints.
    double epsilon = 1e-7;
    bool inX = px >= std::min(u.x, v.x) - epsilon && px <= std::max(u.x, v.x) + epsilon;
    bool inY = py >= std::min(u.y, v.y) - epsilon && py <= std::max(u.y, v.y) + epsilon;
    
    return inX && inY;
}