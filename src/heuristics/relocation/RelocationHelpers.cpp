#include "RelocationHelpers.hpp"

#include <algorithm>
#include <cmath>

namespace {
const int INSIDE = 0;
const int LEFT = 1;
const int RIGHT = 2;
const int BOTTOM = 4;
const int TOP = 8;
} // namespace

namespace relocation_detail {

int computeOutCode(double x, double y, const RegionOfInterest& roi) {
    int code = INSIDE;
    if (x < roi.minX) code |= LEFT;
    else if (x > roi.maxX) code |= RIGHT;

    if (y < roi.minY) code |= BOTTOM;
    else if (y > roi.maxY) code |= TOP;

    return code;
}

int toColIndex(double x, const SpatialGrid& grid) {
    if (grid.getNumCellsX() <= 0) return 0;
    int col = static_cast<int>(std::floor((x - grid.getMinX()) / grid.getCellWidth()));
    return std::clamp(col, 0, grid.getNumCellsX() - 1);
}

int toRowIndex(double y, const SpatialGrid& grid) {
    if (grid.getNumCellsY() <= 0) return 0;
    int row = static_cast<int>(std::floor((y - grid.getMinY()) / grid.getCellHeight()));
    return std::clamp(row, 0, grid.getNumCellsY() - 1);
}

int toCellIndex(int col, int row, const SpatialGrid& grid) {
    return row * grid.getNumCellsX() + col;
}

std::vector<int> collectOriginalNeighbors(int nodeId, const PlanarizedGraph& pGraph) {
    std::vector<int> originalNeighbors;
    std::unordered_set<int> seenNeighbors;
    if (!pGraph.hasNode(nodeId)) return originalNeighbors;

    const auto& startNode = pGraph.getNode(nodeId);

    for (int neighborId : startNode.adjacent_planar_nodes) {
        int current = neighborId;
        int previous = nodeId;

        while (true) {
            if (!pGraph.hasNode(current)) break;

            const auto& currNode = pGraph.getNode(current);

            if (currNode.type == PlanarizedGraph::NodeType::ORIGINAL) {
                if (seenNeighbors.insert(current).second) {
                    originalNeighbors.push_back(current);
                }
                break;
            }

            int next = -1;
            if (currNode.e1_neighbor_prev == previous) next = currNode.e1_neighbor_next;
            else if (currNode.e1_neighbor_next == previous) next = currNode.e1_neighbor_prev;
            else if (currNode.e2_neighbor_prev == previous) next = currNode.e2_neighbor_next;
            else if (currNode.e2_neighbor_next == previous) next = currNode.e2_neighbor_prev;

            if (next == -1) break;

            previous = current;
            current = next;
        }
    }

    return originalNeighbors;
}

std::unordered_set<int> collectIncidentOriginalEdgeIds(int nodeId, const PlanarizedGraph& pGraph) {
    std::unordered_set<int> incidentOriginalEdgeIds;

    if (!pGraph.hasNode(nodeId)) return incidentOriginalEdgeIds;

    // If you have adjacency, use it here instead of forEachEdge over all edges.
    pGraph.forEachEdge([&](int /*edgeId*/, const PlanarizedGraph::PlanarEdge& planarEdge) {
        if (planarEdge.u_id == nodeId || planarEdge.v_id == nodeId) {
            incidentOriginalEdgeIds.insert(planarEdge.original_edge_id);
        }
    });

    return incidentOriginalEdgeIds;
}

bool pointInPolygon(double x, double y, const std::vector<FaceVertex>& polygon) {
    if (polygon.size() < 3) return false;

    bool inside = false;
    for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const auto& pi = polygon[i];
        const auto& pj = polygon[j];

        const bool intersects = ((pi.y > y) != (pj.y > y)) &&
                                (x < (pj.x - pi.x) * (y - pi.y) / (pj.y - pi.y) + pi.x);
        if (intersects) inside = !inside;
    }

    return inside;
}

} // namespace relocation_detail