#include "../RelocationManager.hpp"
#include "RelocationHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <random>
#include <unordered_set>

namespace {
// Numerical tolerance for floating-point comparisons
constexpr double EPS = 1e-12;

// Cohen-Sutherland outcode bit constants
constexpr int OUT_INSIDE = 0;
constexpr int OUT_LEFT = 1;
constexpr int OUT_RIGHT = 2;
constexpr int OUT_BOTTOM = 4;
constexpr int OUT_TOP = 8;
} // namespace

RelocationManager::RelocationManager(PlanarizedGraph& pGraph)
    : m_pGraph(pGraph) {
    std::random_device rd;
    m_rng.seed(rd());  // Seed RNG once in constructor
    
    m_pGraph.forEachNode([&](int, const PlanarizedGraph::PlanarNode& node) {
        if (node.type == PlanarizedGraph::NodeType::ORIGINAL) {
            m_originalNodeIds.push_back(node.id);
            m_originalNodeIdSet.insert(node.id);  // Populate membership set
        }
    });
}

int RelocationManager::selectVariableNode() {
    if (m_originalNodeIds.empty()) return -1;

    if (const char* forced = std::getenv("COPILOT_DEBUG_FORCE_NODE")) {
        char* end = nullptr;
        const long parsed = std::strtol(forced, &end, 10);
        if (end != forced && *end == '\0') {
            const int forcedId = static_cast<int>(parsed);
            if (m_originalNodeIdSet.count(forcedId) > 0) {  // O(1) membership test
                return forcedId;
            }
        }
    }

    if (const char* legacy = std::getenv("COPILOT_DEBUG_FORCE_NODE4")) {
        if (legacy[0] != '\0' && m_originalNodeIdSet.count(4) > 0) {  // O(1) membership test
            return 4;
        }
    }

    std::uniform_int_distribution<> distrib(0, m_originalNodeIds.size() - 1);
    return m_originalNodeIds[distrib(m_rng)];  // Use seeded member RNG
}

RegionOfInterest RelocationManager::calculateROI(int nodeId) {
    RegionOfInterest roi;
    const auto& grid = m_pGraph.getGrid();

    if (!m_pGraph.hasNode(nodeId)) {
        roi.minCol = roi.maxCol = 0;
        roi.minRow = roi.maxRow = 0;
        roi.minX = roi.maxX = grid.getMinX();
        roi.minY = roi.maxY = grid.getMinY();
        return roi;
    }

    const auto& node = m_pGraph.getNode(nodeId);
    auto neighbors = relocation_detail::collectOriginalNeighbors(nodeId, m_pGraph);

    roi.minCol = relocation_detail::toColIndex(node.x, grid);
    roi.maxCol = roi.minCol;
    roi.minRow = relocation_detail::toRowIndex(node.y, grid);
    roi.maxRow = roi.minRow;

    for (int neighborId : neighbors) {
        if (!m_pGraph.hasNode(neighborId)) continue;

        const auto& neighbor = m_pGraph.getNode(neighborId);
        int col = relocation_detail::toColIndex(neighbor.x, grid);
        int row = relocation_detail::toRowIndex(neighbor.y, grid);

        roi.minCol = std::min(roi.minCol, col);
        roi.maxCol = std::max(roi.maxCol, col);
        roi.minRow = std::min(roi.minRow, row);
        roi.maxRow = std::max(roi.maxRow, row);
    }

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
    const auto incidentOriginalEdgeIds = relocation_detail::collectIncidentOriginalEdgeIds(variableNodeId, m_pGraph);

    // Estimate and reserve capacity: 4 border segments + cells * avg edges per cell + rays
    int numCells = (roi.maxCol - roi.minCol + 1) * (roi.maxRow - roi.minRow + 1);
    int estimatedCapacity = 4 + (numCells * 8) + 32;  // 4 borders + edges + rays
    localGeometry.reserve(estimatedCapacity);
    processedEdgeIds.reserve(numCells * 8);

    localGeometry.push_back({roi.minX, roi.minY, roi.maxX, roi.minY, -1, LocalSegmentType::ROI_BORDER});
    localGeometry.push_back({roi.minX, roi.maxY, roi.maxX, roi.maxY, -1, LocalSegmentType::ROI_BORDER});
    localGeometry.push_back({roi.minX, roi.minY, roi.minX, roi.maxY, -1, LocalSegmentType::ROI_BORDER});
    localGeometry.push_back({roi.maxX, roi.minY, roi.maxX, roi.maxY, -1, LocalSegmentType::ROI_BORDER});

    for (int col = roi.minCol; col <= roi.maxCol; ++col) {
        for (int row = roi.minRow; row <= roi.maxRow; ++row) {
            int cellIndex = relocation_detail::toCellIndex(col, row, grid);
            const auto& edgeIdsInCell = grid.getCell(cellIndex).edgeIndices;

            for (int edgeId : edgeIdsInCell) {
                if (!processedEdgeIds.insert(edgeId).second) continue;

                if (!m_pGraph.hasEdge(edgeId)) continue;
                const auto& edge = m_pGraph.getEdge(edgeId);

                if (incidentOriginalEdgeIds.count(edge.original_edge_id) > 0) continue;

                if (!m_pGraph.hasNode(edge.u_id) || !m_pGraph.hasNode(edge.v_id)) continue;

                const auto& u = m_pGraph.getNode(edge.u_id);
                const auto& v = m_pGraph.getNode(edge.v_id);

                auto clippedSegment = clipToBoundingBox(u.x, u.y, v.x, v.y, roi, edge.original_edge_id);
                if (clippedSegment.has_value()) {
                    localGeometry.emplace_back(std::move(clippedSegment.value()));  // Move to avoid copy
                }
            }
        }
    }

    appendRaysToLocalGeometry(localGeometry, roi, variableNodeId);

    return localGeometry;
}

std::optional<LocalSegment> RelocationManager::clipToBoundingBox(
    double x1, double y1, double x2, double y2, const RegionOfInterest& roi, int edgeId) {

    int outcode0 = relocation_detail::computeOutCode(x1, y1, roi);
    int outcode1 = relocation_detail::computeOutCode(x2, y2, roi);
    bool accept = false;

    while (true) {
        if (!(outcode0 | outcode1)) {
            accept = true;
            break;
        } else if (outcode0 & outcode1) {
            break;
        } else {
            double x, y;
            int outcodeOut = outcode0 ? outcode0 : outcode1;

            // Handle clipping to each boundary with explicit parallel segment checks
            if (outcodeOut & OUT_TOP) {
                double dy = y2 - y1;
                if (std::abs(dy) > EPS) {
                    // Segment not parallel to horizontal boundary: compute intersection
                    x = x1 + (x2 - x1) * (roi.maxY - y1) / dy;
                } else {
                    // Segment parallel to horizontal boundary; use appropriate endpoint
                    x = (outcodeOut == outcode0) ? x1 : x2;
                }
                y = roi.maxY;
            } else if (outcodeOut & OUT_BOTTOM) {
                double dy = y2 - y1;
                if (std::abs(dy) > EPS) {
                    // Segment not parallel to horizontal boundary: compute intersection
                    x = x1 + (x2 - x1) * (roi.minY - y1) / dy;
                } else {
                    // Segment parallel to horizontal boundary; use appropriate endpoint
                    x = (outcodeOut == outcode0) ? x1 : x2;
                }
                y = roi.minY;
            } else if (outcodeOut & OUT_RIGHT) {
                double dx = x2 - x1;
                if (std::abs(dx) > EPS) {
                    // Segment not parallel to vertical boundary: compute intersection
                    y = y1 + (y2 - y1) * (roi.maxX - x1) / dx;
                } else {
                    // Segment parallel to vertical boundary; use appropriate endpoint
                    y = (outcodeOut == outcode0) ? y1 : y2;
                }
                x = roi.maxX;
            } else if (outcodeOut & OUT_LEFT) {
                double dx = x2 - x1;
                if (std::abs(dx) > EPS) {
                    // Segment not parallel to vertical boundary: compute intersection
                    y = y1 + (y2 - y1) * (roi.minX - x1) / dx;
                } else {
                    // Segment parallel to vertical boundary; use appropriate endpoint
                    y = (outcodeOut == outcode0) ? y1 : y2;
                }
                x = roi.minX;
            }

            if (outcodeOut == outcode0) {
                x1 = x;
                y1 = y;
                outcode0 = relocation_detail::computeOutCode(x1, y1, roi);
            } else {
                x2 = x;
                y2 = y;
                outcode1 = relocation_detail::computeOutCode(x2, y2, roi);
            }
        }
    }

    if (accept) return LocalSegment{x1, y1, x2, y2, edgeId, LocalSegmentType::ORIGINAL_SUBSEGMENT};
    return std::nullopt;
}

void RelocationManager::appendRaysToLocalGeometry(std::vector<LocalSegment>& localGeometry,
                                                  const RegionOfInterest& roi,
                                                  int variableNodeId) {
    auto neighbors = relocation_detail::collectOriginalNeighbors(variableNodeId, m_pGraph);

    m_pGraph.forEachNode([&](int vId, const PlanarizedGraph::PlanarNode& vNode) {
        if (vNode.type != PlanarizedGraph::NodeType::ORIGINAL) return;
        if (vId == variableNodeId) return;

        if (vNode.x < roi.minX || vNode.x > roi.maxX ||
            vNode.y < roi.minY || vNode.y > roi.maxY) {
            return;
        }

        for (int neighborId : neighbors) {
            if (vId == neighborId) continue;

            if (!m_pGraph.hasNode(neighborId)) continue;
            const auto& xNode = m_pGraph.getNode(neighborId);

            double dx = vNode.x - xNode.x;
            double dy = vNode.y - xNode.y;

            if (std::abs(dx) < EPS && std::abs(dy) < EPS) continue;

            double tMin = std::numeric_limits<double>::infinity();

            if (dx > EPS) tMin = std::min(tMin, (roi.maxX - vNode.x) / dx);
            else if (dx < -EPS) tMin = std::min(tMin, (roi.minX - vNode.x) / dx);

            if (dy > EPS) tMin = std::min(tMin, (roi.maxY - vNode.y) / dy);
            else if (dy < -EPS) tMin = std::min(tMin, (roi.minY - vNode.y) / dy);

            if (std::isinf(tMin) || tMin <= EPS) continue;

            double endX = vNode.x + tMin * dx;
            double endY = vNode.y + tMin * dy;

            LocalSegment ray;
            ray.x1 = vNode.x;
            ray.y1 = vNode.y;
            ray.x2 = endX;
            ray.y2 = endY;
            ray.originalEdgeId = -1;
            ray.type = LocalSegmentType::RAY;
            ray.raySourceNodeId = vId;
            localGeometry.push_back(ray);
        }
    });
}
