#include "RelocationManager.hpp"
#include <random>
#include <unordered_set>
#include <algorithm>
#include <cmath>

// Anonymous namespace for internal helper functions
namespace {
    // Cohen-Sutherland Outcodes
    const int INSIDE = 0; // 0000
    const int LEFT   = 1; // 0001
    const int RIGHT  = 2; // 0010
    const int BOTTOM = 4; // 0100
    const int TOP    = 8; // 1000

    // Helper to compute the region code for a point (x, y)
    int computeOutCode(double x, double y, const RegionOfInterest& roi) {
        int code = INSIDE;
        if (x < roi.minX)      code |= LEFT;
        else if (x > roi.maxX) code |= RIGHT;
        
        if (y < roi.minY)      code |= BOTTOM;
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

    /**
     * @brief Follows planar edges through crossing nodes to find the endpoint ORIGINAL nodes.
     */
    std::vector<int> collectOriginalNeighbors(int nodeId, const PlanarizedGraph& pGraph) {
        std::vector<int> originalNeighbors;
        auto it = pGraph.nodes.find(nodeId);
        if (it == pGraph.nodes.end()) return originalNeighbors;

        const auto& startNode = it->second;

        // Since the variable node is always ORIGINAL, we iterate its incident planar edges
        for (int neighborId : startNode.adjacent_planar_nodes) {
            int current = neighborId;
            int previous = nodeId;

            // Walk the path of the edge
            while (true) {
                auto currIt = pGraph.nodes.find(current);
                if (currIt == pGraph.nodes.end()) break;

                const auto& currNode = currIt->second;
                
                // If we found an ORIGINAL node, we've found the true neighbor
                if (currNode.type == PlanarizedGraph::NodeType::ORIGINAL) {
                    originalNeighbors.push_back(current);
                    break;
                }

                // If it's a CROSSING, we must follow the edge it belongs to.
                // A crossing has two edges (e1 and e2). We check which one we came from.
                int next = -1;
                if (currNode.e1_neighbor_prev == previous)      next = currNode.e1_neighbor_next;
                else if (currNode.e1_neighbor_next == previous) next = currNode.e1_neighbor_prev;
                else if (currNode.e2_neighbor_prev == previous) next = currNode.e2_neighbor_next;
                else if (currNode.e2_neighbor_next == previous) next = currNode.e2_neighbor_prev;

                if (next == -1) break; // Should not happen in a valid planarized structure

                previous = current;
                current = next;
            }
        }
        return originalNeighbors;
    }
}

RelocationManager::RelocationManager(PlanarizedGraph& pGraph)
    : m_pGraph(pGraph) {
    for (const auto& pair : m_pGraph.nodes) {
        const auto& node = pair.second;
        if (node.type == PlanarizedGraph::NodeType::ORIGINAL) {
            m_originalNodeIds.push_back(node.id);
        }
    }
}

int RelocationManager::selectVariableNode() {
    if (m_originalNodeIds.empty()) return -1;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, m_originalNodeIds.size() - 1);

    return m_originalNodeIds[distrib(gen)];
}

RegionOfInterest RelocationManager::calculateROI(int nodeId) {
    RegionOfInterest roi;
    const auto& grid = m_pGraph.getGrid();

    auto nodeIt = m_pGraph.nodes.find(nodeId);
    if (nodeIt == m_pGraph.nodes.end()) {
        roi.minCol = roi.maxCol = 0;
        roi.minRow = roi.maxRow = 0;
        roi.minX = roi.maxX = grid.getMinX();
        roi.minY = roi.maxY = grid.getMinY();
        return roi;
    }
    
    const auto& node = nodeIt->second;
    
    // 1. Get the TRUE ORIGINAL neighbors by traversing through crossings
    auto neighbors = collectOriginalNeighbors(nodeId, m_pGraph);
    
    // 2. Initialize with the variable node's own grid position
    roi.minCol = toColIndex(node.x, grid);
    roi.maxCol = roi.minCol;
    roi.minRow = toRowIndex(node.y, grid);
    roi.maxRow = roi.minRow;

    // 3. Expand grid indices to cover all original neighbors
    for (int neighborId : neighbors) {
        auto neighborIt = m_pGraph.nodes.find(neighborId);
        if (neighborIt == m_pGraph.nodes.end()) continue;

        const auto& neighbor = neighborIt->second;
        int col = toColIndex(neighbor.x, grid);
        int row = toRowIndex(neighbor.y, grid);

        roi.minCol = std::min(roi.minCol, col);
        roi.maxCol = std::max(roi.maxCol, col);
        roi.minRow = std::min(roi.minRow, row);
        roi.maxRow = std::max(roi.maxRow, row);
    }

    // 4. Snap the physical bounds to the outer edges of these grid cells
    roi.minX = grid.getMinX() + (roi.minCol * grid.getCellWidth());
    roi.maxX = grid.getMinX() + ((roi.maxCol + 1) * grid.getCellWidth());
    roi.minY = grid.getMinY() + (roi.minRow * grid.getCellHeight());
    roi.maxY = grid.getMinY() + ((roi.maxRow + 1) * grid.getCellHeight());

    return roi;
}

std::vector<LocalSegment> RelocationManager::extractAndClipGeometry(const RegionOfInterest& roi, int variableNodeId) {
    std::vector<LocalSegment> localGeometry;
    std::unordered_set<int> processedEdgeIds;
    const auto& grid = m_pGraph.getGrid();

    // ROI borders as "Static Walls"
    localGeometry.push_back({roi.minX, roi.minY, roi.maxX, roi.minY, -1}); // Bottom
    localGeometry.push_back({roi.minX, roi.maxY, roi.maxX, roi.maxY, -1}); // Top
    localGeometry.push_back({roi.minX, roi.minY, roi.minX, roi.maxY, -1}); // Left
    localGeometry.push_back({roi.maxX, roi.minY, roi.maxX, roi.maxY, -1}); // Right

    for (int col = roi.minCol; col <= roi.maxCol; ++col) {
        for (int row = roi.minRow; row <= roi.maxRow; ++row) {
            int cellIndex = toCellIndex(col, row, grid);
            const auto& edgeIdsInCell = grid.getCell(cellIndex).edgeIndices;

            for (int edgeId : edgeIdsInCell) {
                if (!processedEdgeIds.insert(edgeId).second) continue; 

                auto edgeIt = m_pGraph.edges.find(edgeId);
                if (edgeIt == m_pGraph.edges.end()) continue;
                const auto& edge = edgeIt->second;

                // Ignore edges connected to the moving node
                if (edge.u_id == variableNodeId || edge.v_id == variableNodeId) continue;

                auto uIt = m_pGraph.nodes.find(edge.u_id);
                auto vIt = m_pGraph.nodes.find(edge.v_id);
                if (uIt == m_pGraph.nodes.end() || vIt == m_pGraph.nodes.end()) continue;

                const auto& u = uIt->second;
                const auto& v = vIt->second;

                auto clippedSegment = clipToBoundingBox(u.x, u.y, v.x, v.y, roi, edgeId);
                if (clippedSegment.has_value()) {
                    localGeometry.push_back(clippedSegment.value());
                }
            }
        }
    }
    return localGeometry;
}

std::optional<LocalSegment> RelocationManager::clipToBoundingBox(
    double x1, double y1, double x2, double y2, const RegionOfInterest& roi, int edgeId) {
    
    int outcode0 = computeOutCode(x1, y1, roi);
    int outcode1 = computeOutCode(x2, y2, roi);
    bool accept = false;

    while (true) {
        if (!(outcode0 | outcode1)) { accept = true; break; }
        else if (outcode0 & outcode1) { break; }
        else {
            double x, y;
            int outcodeOut = outcode1 > outcode0 ? outcode1 : outcode0;

            if (outcodeOut & TOP) {
                x = x1 + (x2 - x1) * (roi.maxY - y1) / (y2 - y1);
                y = roi.maxY;
            } else if (outcodeOut & BOTTOM) {
                x = x1 + (x2 - x1) * (roi.minY - y1) / (y2 - y1);
                y = roi.minY;
            } else if (outcodeOut & RIGHT) {
                y = y1 + (y2 - y1) * (roi.maxX - x1) / (x2 - x1);
                x = roi.maxX;
            } else if (outcodeOut & LEFT) {
                y = y1 + (y2 - y1) * (roi.minX - x1) / (x2 - x1);
                x = roi.minX;
            }

            if (outcodeOut == outcode0) {
                x1 = x; y1 = y;
                outcode0 = computeOutCode(x1, y1, roi);
            } else {
                x2 = x; y2 = y;
                outcode1 = computeOutCode(x2, y2, roi);
            }
        }
    }

    if (accept) return LocalSegment{x1, y1, x2, y2, edgeId};
    return std::nullopt;
}