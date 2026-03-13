#include "intersection_detector.hpp"
#include <unordered_set>
#include <cmath>
#include <algorithm>

// Helper to create a unique 64-bit key for any pair of 32-bit edge IDs.
// This ensures we don't check Edge A vs Edge B multiple times if they share multiple cells.
inline uint64_t makeEdgePairKey(int id1, int id2) {
    uint64_t min_id = static_cast<uint64_t>(std::min(id1, id2));
    uint64_t max_id = static_cast<uint64_t>(std::max(id1, id2));
    return (min_id << 32) | max_id;
}

std::vector<PlanarizedGraph::IntersectionData> findIntersections(const Graph& graph, const SpatialGrid& grid) {
    std::vector<PlanarizedGraph::IntersectionData> intersections;
    std::unordered_set<uint64_t> checkedPairs;
    
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
                const auto& e1 = graph.edges[edge1_id];
                const auto& e2 = graph.edges[edge2_id];

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