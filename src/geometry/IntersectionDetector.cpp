#include "IntersectionDetector.hpp"
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <cstdint>

// Helper to create a unique 64-bit key for any pair of 32-bit edge IDs.
// This ensures we don't check Edge A vs Edge B multiple times if they share multiple cells.
inline uint64_t makeEdgePairKey(int id1, int id2) {
    uint64_t min_id = static_cast<uint64_t>(std::min(id1, id2));
    uint64_t max_id = static_cast<uint64_t>(std::max(id1, id2));
    return (min_id << 32) | max_id;
}

std::vector<PlanarizedGraph::IntersectionData> IntersectionDetector::findIntersections(const Graph& graph, const SpatialGrid& grid) {
    std::vector<PlanarizedGraph::IntersectionData> intersections;
    std::unordered_set<uint64_t> checkedPairs;

    // Graph::Edge IDs are not guaranteed to be contiguous vector indices.
    // Build an ID -> edge lookup to avoid out-of-bounds vector access.
    std::unordered_map<int, const Graph::Edge*> edgeById;
    edgeById.reserve(graph.edges.size());
    for (const auto& edge : graph.edges) {
        edgeById[edge.id] = &edge;
    }
    
    // We use a small epsilon to avoid floating point errors identifying 
    // shared endpoints (nodes) as crossings.
    const double EPS = 1e-9; 

    // 1. Iterate over every cell in the spatial grid
    size_t numCells = grid.getNumCells();
    for (size_t c = 0; c < numCells; ++c) {
        const auto& cell = grid.getCell(c);
        const auto& cellEdges = cell.edgeIndices;
        size_t numEdgesInCell = cellEdges.size();

        // 2. Pairwise comparison of edges within this specific cell
        for (size_t i = 0; i < numEdgesInCell; ++i) {
            for (size_t j = i + 1; j < numEdgesInCell; ++j) {
                int edge1_id = cellEdges[i];
                int edge2_id = cellEdges[j];

                // Skip if we've already checked this pair in a different cell
                uint64_t pairKey = makeEdgePairKey(edge1_id, edge2_id);
                if (!checkedPairs.insert(pairKey).second) {
                    continue; // .second is false if the key was already in the set
                }

                // 3. Fetch the actual graph edges and their geometric nodes
                auto e1It = edgeById.find(edge1_id);
                auto e2It = edgeById.find(edge2_id);
                if (e1It == edgeById.end() || e2It == edgeById.end()) {
                    continue;
                }

                const auto& e1 = *(e1It->second);
                const auto& e2 = *(e2It->second);

                // Optimization: If they share an original node, they can't cross in the middle.
                if (e1.u_id == e2.u_id || e1.u_id == e2.v_id || 
                    e1.v_id == e2.u_id || e1.v_id == e2.v_id) {
                    continue;
                }

                const auto& A = graph.nodes[e1.u_id];
                const auto& B = graph.nodes[e1.v_id];
                const auto& C = graph.nodes[e2.u_id];
                const auto& D = graph.nodes[e2.v_id];

                // 4. Vector Cross-Product Math
                // r = B - A (Direction vector of Edge 1)
                double rx = B.x - A.x;
                double ry = B.y - A.y;
                
                // s = D - C (Direction vector of Edge 2)
                double sx = D.x - C.x;
                double sy = D.y - C.y;

                // Cross product: r x s
                double denom = (rx * sy) - (ry * sx);

                // If denominator is 0, lines are parallel or collinear
                if (std::abs(denom) < EPS) {
                    continue; 
                }

                // diff = C - A
                double diff_x = C.x - A.x;
                double diff_y = C.y - A.y;

                // t1 = (diff x s) / denom
                double t1 = ((diff_x * sy) - (diff_y * sx)) / denom;
                
                // t2 = (diff x r) / denom
                double t2 = ((diff_x * ry) - (diff_y * rx)) / denom;

                // 5. Check if the intersection happens STRICTLY inside both segments
                if (t1 > EPS && t1 < (1.0 - EPS) && t2 > EPS && t2 < (1.0 - EPS)) {
                    // Calculate exact X and Y
                    double cross_x = A.x + (t1 * rx);
                    double cross_y = A.y + (t1 * ry);

                    intersections.push_back({
                        cross_x, 
                        cross_y, 
                        edge1_id, 
                        edge2_id, 
                        t1, 
                        t2
                    });
                }
            }
        }
    }

    return intersections;
}