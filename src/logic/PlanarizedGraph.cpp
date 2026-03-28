#include "PlanarizedGraph.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

PlanarizedGraph::PlanarizedGraph(const Graph& originalGraph, const std::vector<IntersectionData>& intersections) {
    
    // 1. Calculate boundaries to initialize the grid
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& node : originalGraph.nodes) {
        if (node.x < minX) minX = node.x;
        if (node.x > maxX) maxX = node.x;
        if (node.y < minY) minY = node.y;
        if (node.y > maxY) maxY = node.y;
    }

    // Initialize the member grid with actual dimensions
    int totalExpectedNodes = originalGraph.nodes.size() + intersections.size();
    grid = SpatialGrid(minX, maxX, minY, maxY, totalExpectedNodes);

    // 2. Map ORIGINAL nodes
    int maxOrigNodeId = -1;
    for (const auto& origNode : originalGraph.nodes) {
        PlanarNode pNode;
        pNode.id = origNode.id; pNode.x = origNode.x; pNode.y = origNode.y;
        pNode.type = NodeType::ORIGINAL;
        pNode.original_node_id = origNode.id;
        
        nodes[pNode.id] = pNode;
        grid.insertNode(pNode.id, pNode.x, pNode.y); // <-- Sync to grid
        
        if (origNode.id > maxOrigNodeId) maxOrigNodeId = origNode.id;
    }
    
    nextNodeId = maxOrigNodeId + 1;
    
    // 3. Setup Sort-by-t Buckets
    std::unordered_map<int, std::vector<std::pair<double, int>>> edgeBuckets;
    for (const auto& origEdge : originalGraph.edges) {
        edgeBuckets[origEdge.id].push_back({0.0, origEdge.u_id});
        edgeBuckets[origEdge.id].push_back({1.0, origEdge.v_id});
    }
    
    // 4. Create CROSSING nodes
    for (const auto& ix : intersections) {
        int cId = getNextNodeId();
        
        PlanarNode cNode;
        cNode.id = cId; cNode.x = ix.x; cNode.y = ix.y;
        cNode.type = NodeType::CROSSING;
        cNode.original_edge_1 = ix.edge1_id;
        cNode.original_edge_2 = ix.edge2_id;
        
        nodes[cId] = cNode;
        grid.insertNode(cId, cNode.x, cNode.y); // <-- Sync to grid
        
        edgeBuckets[ix.edge1_id].push_back({ix.t1, cId});
        edgeBuckets[ix.edge2_id].push_back({ix.t2, cId});
    }
    
    // 5. Stitch Edges (addPlanarEdge will handle grid insertion automatically)
    for (auto& pair : edgeBuckets) {
        int origEdgeId = pair.first;
        auto& bucket = pair.second;
        
        std::sort(bucket.begin(), bucket.end());
        
        for (size_t i = 0; i < bucket.size() - 1; ++i) {
            int nodeA_id = bucket[i].second;
            int nodeB_id = bucket[i+1].second;
            
            addPlanarEdge(nodeA_id, nodeB_id, origEdgeId); // Grid updated inside here!
            
            // ... (Your exact same logic for linking Node A and Node B pointers) ...
            PlanarNode& nodeA = nodes[nodeA_id];
            if (nodeA.type == NodeType::ORIGINAL) nodeA.adjacent_planar_nodes.push_back(nodeB_id);
            else {
                if (nodeA.original_edge_1 == origEdgeId) nodeA.e1_neighbor_next = nodeB_id;
                else if (nodeA.original_edge_2 == origEdgeId) nodeA.e2_neighbor_next = nodeB_id;
            }
            
            PlanarNode& nodeB = nodes[nodeB_id];
            if (nodeB.type == NodeType::ORIGINAL) nodeB.adjacent_planar_nodes.push_back(nodeA_id);
            else {
                if (nodeB.original_edge_1 == origEdgeId) nodeB.e1_neighbor_prev = nodeA_id;
                else if (nodeB.original_edge_2 == origEdgeId) nodeB.e2_neighbor_prev = nodeA_id;
            }
        }
    }
}

// --- SYNCHRONIZED HELPER METHODS ---

void PlanarizedGraph::addPlanarEdge(int u_id, int v_id, int origEdgeId) {
    int newId = getNextEdgeId();
    edges[newId] = {newId, u_id, v_id, origEdgeId};
    originalEdgeToPlanarEdges[origEdgeId].push_back(newId);

    // Automatically sync the new sub-segment to the spatial grid
    PlanarNode& u = nodes[u_id];
    PlanarNode& v = nodes[v_id];
    grid.insertEdge(newId, u.x, u.y, v.x, v.y); 
}

void PlanarizedGraph::removePlanarEdge(int u_id, int v_id, int origEdgeId) {
    auto& subEdges = originalEdgeToPlanarEdges[origEdgeId];
    for (auto it = subEdges.begin(); it != subEdges.end(); ++it) {
        PlanarEdge& pe = edges[*it];
        if ((pe.u_id == u_id && pe.v_id == v_id) || (pe.u_id == v_id && pe.v_id == u_id)) {
            int edgeIdToDelete = *it;
            
            // Automatically remove the sub-segment from the spatial grid
            PlanarNode& u = nodes[pe.u_id];
            PlanarNode& v = nodes[pe.v_id];
            grid.removeEdge(edgeIdToDelete, u.x, u.y, v.x, v.y);

            edges.erase(edgeIdToDelete);
            subEdges.erase(it);
            break; 
        }
    }
}

// ------------------------------------

void PlanarizedGraph::destroyCrossing(int crossingNodeIdx) {
    PlanarNode& cNode = nodes[crossingNodeIdx];

    // (Your exact healing logic...)
    int e1 = cNode.original_edge_1;
    int prev1 = cNode.e1_neighbor_prev;
    int next1 = cNode.e1_neighbor_next;
    updateNodeNeighbor(prev1, crossingNodeIdx, next1, e1);
    updateNodeNeighbor(next1, crossingNodeIdx, prev1, e1);
    removePlanarEdge(prev1, crossingNodeIdx, e1); // Grid syncs automatically!
    removePlanarEdge(crossingNodeIdx, next1, e1);
    addPlanarEdge(prev1, next1, e1);              // Grid syncs automatically!

    int e2 = cNode.original_edge_2;
    int prev2 = cNode.e2_neighbor_prev;
    int next2 = cNode.e2_neighbor_next;
    updateNodeNeighbor(prev2, crossingNodeIdx, next2, e2);
    updateNodeNeighbor(next2, crossingNodeIdx, prev2, e2);
    removePlanarEdge(prev2, crossingNodeIdx, e2);
    removePlanarEdge(crossingNodeIdx, next2, e2);
    addPlanarEdge(prev2, next2, e2);

    // Sync node removal to grid, then delete from graph
    grid.removeNode(crossingNodeIdx, cNode.x, cNode.y); 
    nodes.erase(crossingNodeIdx); 
}

void PlanarizedGraph::createCrossing(int edge1Id, int edge2Id, double x, double y) {
    int cId = getNextNodeId();
    
    PlanarNode cNode;
    cNode.id = cId; cNode.x = x; cNode.y = y;
    cNode.type = NodeType::CROSSING;
    cNode.original_edge_1 = edge1Id;
    cNode.original_edge_2 = edge2Id;

    int affectedEdges[2] = {edge1Id, edge2Id};

    for (int eId : affectedEdges) {
        int targetEdgeId = -1;
        int u = -1, v = -1;

        for (int pEdgeId : originalEdgeToPlanarEdges[eId]) {
            PlanarEdge& pe = edges[pEdgeId];
            if (isPointOnSegment(x, y, pe.u_id, pe.v_id)) {
                targetEdgeId = pEdgeId;
                u = pe.u_id; v = pe.v_id;
                break; 
            }
        }
        if (targetEdgeId == -1) continue; 

        if (eId == edge1Id) {
            cNode.e1_neighbor_prev = u; cNode.e1_neighbor_next = v;
        } else {
            cNode.e2_neighbor_prev = u; cNode.e2_neighbor_next = v;
        }

        updateNodeNeighbor(u, v, cId, eId);
        updateNodeNeighbor(v, u, cId, eId);

        removePlanarEdge(u, v, eId);   // Grid syncs automatically!
        addPlanarEdge(u, cId, eId);    // Grid syncs automatically!
        addPlanarEdge(cId, v, eId);    // Grid syncs automatically!
    }

    nodes[cId] = cNode;
    grid.insertNode(cId, cNode.x, cNode.y); // Sync node creation to grid
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