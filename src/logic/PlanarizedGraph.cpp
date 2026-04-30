#include "PlanarizedGraph.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

constexpr double kIntersectionEpsilon = 1e-9;

inline std::uint64_t makeCrossingPairKey(int edgeA, int edgeB) {
    const std::uint64_t a = static_cast<std::uint64_t>(std::min(edgeA, edgeB));
    const std::uint64_t b = static_cast<std::uint64_t>(std::max(edgeA, edgeB));
    return (a << 32) | b;
}

} // namespace

void PlanarizedGraph::ensureNodeCapacity(int id) {
    if (id < 0) return;
    if (id >= static_cast<int>(nodes.size())) {
        nodes.resize(id + 1);
        nodeActive.resize(id + 1, 0);
    }
}

void PlanarizedGraph::ensureEdgeCapacity(int id) {
    if (id < 0) return;
    if (id >= static_cast<int>(edges.size())) {
        edges.resize(id + 1);
        edgeActive.resize(id + 1, 0);
    }
}

void PlanarizedGraph::ensureOriginalEdgeCapacity(int id) {
    if (id < 0) return;
    if (id >= static_cast<int>(originalEdgeToPlanarEdges.size())) {
        originalEdgeToPlanarEdges.resize(id + 1);
    }
}

void PlanarizedGraph::deactivateNode(int id) {
    if (!hasNode(id)) return;
    nodeActive[id] = 0;
    freeNodeIds.push_back(id);
}

void PlanarizedGraph::deactivateEdge(int id) {
    if (!hasEdge(id)) return;
    edgeActive[id] = 0;
    freeEdgeIds.push_back(id);
}

int PlanarizedGraph::getNextNodeId() {
    int id;
    if (!freeNodeIds.empty()) {
        id = freeNodeIds.back();
        freeNodeIds.pop_back();
    } else {
        id = nextNodeId++;
    }

    ensureNodeCapacity(id);
    nodes[id] = PlanarNode{};
    nodes[id].id = id;
    nodeActive[id] = 1;
    return id;
}

int PlanarizedGraph::getNextEdgeId() {
    if (!freeEdgeIds.empty()) {
        const int id = freeEdgeIds.back();
        freeEdgeIds.pop_back();
        ensureEdgeCapacity(id);
        edgeActive[id] = 1;
        return id;
    }

    const int id = nextEdgeId++;
    ensureEdgeCapacity(id);
    edgeActive[id] = 1;
    return id;
}

bool PlanarizedGraph::hasNode(int nodeId) const {
    return nodeId >= 0 && nodeId < static_cast<int>(nodeActive.size()) && nodeActive[nodeId] != 0;
}

bool PlanarizedGraph::hasEdge(int edgeId) const {
    return edgeId >= 0 && edgeId < static_cast<int>(edgeActive.size()) && edgeActive[edgeId] != 0;
}

PlanarizedGraph::PlanarNode& PlanarizedGraph::getNode(int nodeId) {
    return nodes[nodeId];
}

const PlanarizedGraph::PlanarNode& PlanarizedGraph::getNode(int nodeId) const {
    return nodes[nodeId];
}

PlanarizedGraph::PlanarEdge& PlanarizedGraph::getEdge(int edgeId) {
    return edges[edgeId];
}

const PlanarizedGraph::PlanarEdge& PlanarizedGraph::getEdge(int edgeId) const {
    return edges[edgeId];
}

PlanarizedGraph::PlanarizedGraph(const Graph& originalGraph, const std::vector<IntersectionData>& intersections) {
    const int totalExpectedNodes = static_cast<int>(originalGraph.nodes.size() + intersections.size());
    const int totalExpectedEdges = static_cast<int>(originalGraph.edges.size() + (2 * intersections.size()));

    int maxOrigNodeId = -1;

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
        if (node.id > maxOrigNodeId) maxOrigNodeId = node.id;
    }

    const int initialNodeSlots = std::max(totalExpectedNodes, maxOrigNodeId + 1 + static_cast<int>(intersections.size()));
    nodes.resize(initialNodeSlots);
    nodeActive.resize(initialNodeSlots, 0);

    edges.reserve(totalExpectedEdges);
    edgeActive.reserve(totalExpectedEdges);
    crossingPairToNode.reserve(intersections.size() * 2);

    // Initialize the member grid with actual dimensions
    grid = SpatialGrid(minX, maxX, minY, maxY, totalExpectedNodes, totalExpectedEdges);

    // 2. Map ORIGINAL nodes
    for (const auto& origNode : originalGraph.nodes) {
        PlanarNode pNode;
        pNode.id = origNode.id; pNode.x = origNode.x; pNode.y = origNode.y;
        pNode.type = NodeType::ORIGINAL;
        pNode.original_node_id = origNode.id;

        ensureNodeCapacity(pNode.id);
        nodes[pNode.id] = pNode;
        nodeActive[pNode.id] = 1;
        grid.insertNode(pNode.id, pNode.x, pNode.y); // <-- Sync to grid
    }
    
    nextNodeId = maxOrigNodeId + 1;

    int maxOrigEdgeId = -1;
    for (const auto& origEdge : originalGraph.edges) {
        if (origEdge.id > maxOrigEdgeId) maxOrigEdgeId = origEdge.id;
    }
    for (const auto& ix : intersections) {
        if (ix.edge1_id > maxOrigEdgeId) maxOrigEdgeId = ix.edge1_id;
        if (ix.edge2_id > maxOrigEdgeId) maxOrigEdgeId = ix.edge2_id;
    }
    if (maxOrigEdgeId >= 0) {
        ensureOriginalEdgeCapacity(maxOrigEdgeId);
    }
    
    // 3. Setup Sort-by-t Buckets
    std::vector<std::vector<std::pair<double, int>>> edgeBuckets(originalEdgeToPlanarEdges.size());
    for (const auto& origEdge : originalGraph.edges) {
        // origEdge.u_id is the index. We need originalGraph.nodes[index].id
        int actual_U_ID = originalGraph.nodes[origEdge.u_id].id;
        int actual_V_ID = originalGraph.nodes[origEdge.v_id].id;

        edgeBuckets[origEdge.id].push_back({0.0, actual_U_ID});
        edgeBuckets[origEdge.id].push_back({1.0, actual_V_ID});
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
        crossingPairToNode[makeCrossingPairKey(ix.edge1_id, ix.edge2_id)] = cId;
        
        ensureOriginalEdgeCapacity(ix.edge1_id);
        ensureOriginalEdgeCapacity(ix.edge2_id);
        if (edgeBuckets.size() < originalEdgeToPlanarEdges.size()) {
            edgeBuckets.resize(originalEdgeToPlanarEdges.size());
        }

        edgeBuckets[ix.edge1_id].push_back({ix.t1, cId});
        edgeBuckets[ix.edge2_id].push_back({ix.t2, cId});
    }
    
    // 5. Stitch Edges (addPlanarEdge will handle grid insertion automatically)
    for (int origEdgeId = 0; origEdgeId < static_cast<int>(edgeBuckets.size()); ++origEdgeId) {
        auto& bucket = edgeBuckets[origEdgeId];
        if (bucket.size() < 2) continue;
        
        std::sort(bucket.begin(), bucket.end());
        
        for (size_t i = 0; i < bucket.size() - 1; ++i) {
            int nodeA_id = bucket[i].second;
            int nodeB_id = bucket[i+1].second;
            
            addPlanarEdge(nodeA_id, nodeB_id, origEdgeId); // Grid updated inside here!
            
            // ... (Your exact same logic for linking Node A and Node B pointers) ...
            PlanarNode& nodeA = getNode(nodeA_id);
            if (nodeA.type == NodeType::ORIGINAL) nodeA.adjacent_planar_nodes.push_back(nodeB_id);
            else {
                if (nodeA.original_edge_1 == origEdgeId) nodeA.e1_neighbor_next = nodeB_id;
                else if (nodeA.original_edge_2 == origEdgeId) nodeA.e2_neighbor_next = nodeB_id;
            }
            
            PlanarNode& nodeB = getNode(nodeB_id);
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
    ensureOriginalEdgeCapacity(origEdgeId);
    originalEdgeToPlanarEdges[origEdgeId].push_back(newId);
    getNode(u_id).incident_planar_edges.push_back(newId);
    getNode(v_id).incident_planar_edges.push_back(newId);

    // Automatically sync the new sub-segment to the spatial grid
    PlanarNode& u = getNode(u_id);
    PlanarNode& v = getNode(v_id);
    grid.insertEdge(newId, u.x, u.y, v.x, v.y); 
}

void PlanarizedGraph::removePlanarEdge(int u_id, int v_id, int origEdgeId) {
    if (origEdgeId < 0 || origEdgeId >= static_cast<int>(originalEdgeToPlanarEdges.size())) return;

    auto& subEdges = originalEdgeToPlanarEdges[origEdgeId];
    for (auto it = subEdges.begin(); it != subEdges.end(); ++it) {
        if (!hasEdge(*it)) continue;
        PlanarEdge& pe = getEdge(*it);
        if ((pe.u_id == u_id && pe.v_id == v_id) || (pe.u_id == v_id && pe.v_id == u_id)) {
            int edgeIdToDelete = *it;

            auto removeFast = [](std::vector<int>& vec, int val) {
                auto found = std::find(vec.begin(), vec.end(), val);
                if (found != vec.end()) {
                    *found = vec.back();
                    vec.pop_back();
                }
            };

            removeFast(getNode(pe.u_id).incident_planar_edges, edgeIdToDelete);
            removeFast(getNode(pe.v_id).incident_planar_edges, edgeIdToDelete);
            
            // Automatically remove the sub-segment from the spatial grid
            PlanarNode& u = getNode(pe.u_id);
            PlanarNode& v = getNode(pe.v_id);
            grid.removeEdge(edgeIdToDelete, u.x, u.y, v.x, v.y);

            deactivateEdge(edgeIdToDelete);
            *it = subEdges.back();
            subEdges.pop_back();
            break; 
        }
    }
}

// ------------------------------------

void PlanarizedGraph::destroyCrossing(int crossingNodeIdx) {
    if (!hasNode(crossingNodeIdx)) return;
    PlanarNode& cNode = getNode(crossingNodeIdx);
    if (cNode.type != NodeType::CROSSING) return;

    const std::uint64_t pairKey = makeCrossingPairKey(cNode.original_edge_1, cNode.original_edge_2);
    crossingPairToNode.erase(pairKey);

    // (Your exact healing logic...)
    int e1 = cNode.original_edge_1;
    int prev1 = cNode.e1_neighbor_prev;
    int next1 = cNode.e1_neighbor_next;
    if (hasNode(prev1) && hasNode(next1)) {
        updateNodeNeighbor(prev1, crossingNodeIdx, next1, e1);
        updateNodeNeighbor(next1, crossingNodeIdx, prev1, e1);
        removePlanarEdge(prev1, crossingNodeIdx, e1); // Grid syncs automatically!
        removePlanarEdge(crossingNodeIdx, next1, e1);
        addPlanarEdge(prev1, next1, e1);              // Grid syncs automatically!
    }

    int e2 = cNode.original_edge_2;
    int prev2 = cNode.e2_neighbor_prev;
    int next2 = cNode.e2_neighbor_next;
    if (hasNode(prev2) && hasNode(next2)) {
        updateNodeNeighbor(prev2, crossingNodeIdx, next2, e2);
        updateNodeNeighbor(next2, crossingNodeIdx, prev2, e2);
        removePlanarEdge(prev2, crossingNodeIdx, e2);
        removePlanarEdge(crossingNodeIdx, next2, e2);
        addPlanarEdge(prev2, next2, e2);
    }

    // Sync node removal to grid, then delete from graph
    grid.removeNode(crossingNodeIdx, cNode.x, cNode.y); 
    deactivateNode(crossingNodeIdx);
}

void PlanarizedGraph::createCrossing(int edge1Id, int edge2Id, double x, double y) {
    if (getCrossingNodeForPair(edge1Id, edge2Id) != -1) {
        return;
    }

    const std::uint64_t pairKey = makeCrossingPairKey(edge1Id, edge2Id);
    int cId = getNextNodeId();

    PlanarNode& storedNode = nodes[cId];
    storedNode = PlanarNode{};
    storedNode.id = cId;
    storedNode.x = x;
    storedNode.y = y;
    storedNode.type = NodeType::CROSSING;
    storedNode.original_edge_1 = edge1Id;
    storedNode.original_edge_2 = edge2Id;
    grid.insertNode(cId, x, y);
    crossingPairToNode[pairKey] = cId;

    auto rollbackCrossingCreation = [&](const std::string& reason) {
        (void)reason;
        crossingPairToNode.erase(pairKey);
        grid.removeNode(cId, x, y);
        deactivateNode(cId);
    };

    struct EdgeSplitTarget {
        int edgeId;
        int u;
        int v;
    };

    std::vector<EdgeSplitTarget> targets;
    targets.reserve(2);
    int affectedEdges[2] = {edge1Id, edge2Id};

    for (int eId : affectedEdges) {
        if (eId < 0 || eId >= static_cast<int>(originalEdgeToPlanarEdges.size())) {
            rollbackCrossingCreation("Invalid original edge ID " + std::to_string(eId));
            return;
        }

        int u = -1;
        int v = -1;
        bool found = false;

        auto& pEdges = originalEdgeToPlanarEdges[eId];

        for (int pEdgeId : pEdges) {
            if (!hasEdge(pEdgeId)) continue;
            PlanarEdge& pe = getEdge(pEdgeId);

            if (isPointOnSegment(x, y, pe.u_id, pe.v_id)) {
                u = pe.u_id;
                v = pe.v_id;
                found = true;
                break;
            }
        }

        if (!found) {
            rollbackCrossingCreation("Point not found on any planar segment of original edge " + std::to_string(eId));
            return;
        }

        targets.push_back({eId, u, v});
    }

    for (const auto& target : targets) {
        const int eId = target.edgeId;
        const int u = target.u;
        const int v = target.v;

        if (eId == edge1Id) {
            storedNode.e1_neighbor_prev = u;
            storedNode.e1_neighbor_next = v;
        } else {
            storedNode.e2_neighbor_prev = u;
            storedNode.e2_neighbor_next = v;
        }

        updateNodeNeighbor(u, v, cId, eId);
        updateNodeNeighbor(v, u, cId, eId);

        removePlanarEdge(u, v, eId);
        addPlanarEdge(u, cId, eId);
        addPlanarEdge(cId, v, eId);
    }
}


void PlanarizedGraph::updateNodeNeighbor(int nodeId, int oldNeighbor, int newNeighbor, int origEdgeId) {
    if (!hasNode(nodeId)) return;
    PlanarNode& node = getNode(nodeId);

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

bool PlanarizedGraph::isPointOnSegment(double px, double py, int u_id, int v_id) {
    PlanarNode& u = getNode(u_id);
    PlanarNode& v = getNode(v_id);

    const double dx = v.x - u.x;
    const double dy = v.y - u.y;
    const double squaredLength = dx * dx + dy * dy;

    if (squaredLength < kIntersectionEpsilon * kIntersectionEpsilon) {
        const double ddx = px - u.x;
        const double ddy = py - u.y;
        return (ddx * ddx + ddy * ddy) <= (kIntersectionEpsilon * kIntersectionEpsilon);
    }

    // Keep the same absolute tolerance as the spatial grid layer so crossing tests agree.
    const double crossProduct = std::abs((px - u.x) * dy - (py - u.y) * dx);
    const double distanceToLine = crossProduct / std::sqrt(squaredLength);
    if (distanceToLine > kIntersectionEpsilon) {
        return false;
    }

    // Check projection is inside segment extents with the same tolerance.
    const double dotProduct = (px - u.x) * dx + (py - u.y) * dy;
    const double eps = kIntersectionEpsilon * squaredLength;
    if (dotProduct < -eps) return false;
    if (dotProduct > squaredLength + eps) return false;

    return true;
}

const SpatialGrid& PlanarizedGraph::getGrid() const {
    return grid;
}

void PlanarizedGraph::updateNodePosition(int nodeId, double newX, double newY) {
    if (!hasNode(nodeId)) return;

    auto& node = getNode(nodeId);
    const double oldX = node.x;
    const double oldY = node.y;

    const std::vector<int>& incidentEdges = node.incident_planar_edges;

    for (int edgeId : incidentEdges) {
        if (!hasEdge(edgeId)) continue;

        const auto& edge = getEdge(edgeId);
        const auto& uNode = getNode(edge.u_id);
        const auto& vNode = getNode(edge.v_id);

        const double x1 = (edge.u_id == nodeId) ? oldX : uNode.x;
        const double y1 = (edge.u_id == nodeId) ? oldY : uNode.y;
        const double x2 = (edge.v_id == nodeId) ? oldX : vNode.x;
        const double y2 = (edge.v_id == nodeId) ? oldY : vNode.y;

        grid.removeEdge(edgeId, x1, y1, x2, y2);
    }

    grid.removeNode(nodeId, oldX, oldY);
    node.x = newX;
    node.y = newY;
    grid.insertNode(nodeId, newX, newY);

    for (int edgeId : incidentEdges) {
        if (!hasEdge(edgeId)) continue;

        const auto& edge = getEdge(edgeId);
        const auto& uNode = getNode(edge.u_id);
        const auto& vNode = getNode(edge.v_id);
        grid.insertEdge(edgeId, uNode.x, uNode.y, vNode.x, vNode.y);
    }
}

int PlanarizedGraph::countTotalCrossings() const {
    int count = 0;
    forEachNode([&count](int, const PlanarNode& node) {
        if (node.type == NodeType::CROSSING) ++count;
    });
    return count;
}

int PlanarizedGraph::getCrossingNodeForPair(int edgeA, int edgeB) const {
    const auto it = crossingPairToNode.find(makeCrossingPairKey(edgeA, edgeB));
    if (it == crossingPairToNode.end()) {
        return -1;
    }
    return hasNode(it->second) ? it->second : -1;
}