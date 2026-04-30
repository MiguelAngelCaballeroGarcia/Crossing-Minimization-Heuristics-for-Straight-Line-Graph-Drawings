#include "../RelocationManager.hpp"
#include "RelocationHelpers.hpp"
#include "../../geometry/FruchtermanReingold_custom.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <boost/concept_check.hpp>

namespace {

constexpr double kIntersectionEpsilon = 1e-9;

// OPTIMIZATION 1 & 2: DDA Grid Traversal and Vector instead of unordered_set
void appendCandidateEdgesForSegmentDDA(const SpatialGrid& grid, 
                                       double x1, double y1, 
                                       double x2, double y2,
                                       std::vector<int>& outCandidates) {
    if (grid.getNumCells() == 0) return;

    const double minX = grid.getMinX();
    const double minY = grid.getMinY();
    const double cellW = grid.getCellWidth();
    const double cellH = grid.getCellHeight();
    const int numX = grid.getNumCellsX();
    const int numY = grid.getNumCellsY();

    // Ray start and end in grid coordinates
    double startColF = (x1 - minX) / cellW;
    double startRowF = (y1 - minY) / cellH;
    double endColF = (x2 - minX) / cellW;
    double endRowF = (y2 - minY) / cellH;

    int currentX = std::clamp(static_cast<int>(std::floor(startColF)), 0, numX - 1);
    int currentY = std::clamp(static_cast<int>(std::floor(startRowF)), 0, numY - 1);
    int endX = std::clamp(static_cast<int>(std::floor(endColF)), 0, numX - 1);
    int endY = std::clamp(static_cast<int>(std::floor(endRowF)), 0, numY - 1);

    const double dx = x2 - x1;
    const double dy = y2 - y1;

    const int stepX = (dx > 0) ? 1 : ((dx < 0) ? -1 : 0);
    const int stepY = (dy > 0) ? 1 : ((dy < 0) ? -1 : 0);

    const double tDeltaX = (stepX != 0) ? std::abs(cellW / dx) : std::numeric_limits<double>::infinity();
    const double tDeltaY = (stepY != 0) ? std::abs(cellH / dy) : std::numeric_limits<double>::infinity();

    double tMaxX = (stepX > 0) ? ((currentX + 1) * cellW + minX - x1) / dx :
                   (stepX < 0) ? ((currentX * cellW + minX - x1) / dx) : std::numeric_limits<double>::infinity();
                   
    double tMaxY = (stepY > 0) ? ((currentY + 1) * cellH + minY - y1) / dy :
                   (stepY < 0) ? ((currentY * cellH + minY - y1) / dy) : std::numeric_limits<double>::infinity();

    // DDA Traversal
    while (true) {
        const int cellIndex = relocation_detail::toCellIndex(currentX, currentY, grid);
        const auto& cell = grid.getCell(cellIndex);
        outCandidates.insert(outCandidates.end(), cell.edgeIndices.begin(), cell.edgeIndices.end());

        if (currentX == endX && currentY == endY) break;

        if (tMaxX < tMaxY) {
            tMaxX += tDeltaX;
            currentX += stepX;
            if (currentX < 0 || currentX >= numX) break;
        } else {
            tMaxY += tDeltaY;
            currentY += stepY;
            if (currentY < 0 || currentY >= numY) break;
        }
    }
}

void collectCandidateOriginalEdgesForSegment(const PlanarizedGraph& pGraph,
                                             const SpatialGrid& grid,
                                             double x1, double y1, double x2, double y2,
                                             std::vector<int>& planarEdgesBuffer,
                                             std::vector<int>& originalEdgesBuffer,
                                             std::vector<bool>& seenTracker) {
    planarEdgesBuffer.clear();
    originalEdgesBuffer.clear();

    appendCandidateEdgesForSegmentDDA(grid, x1, y1, x2, y2, planarEdgesBuffer);

    // Ensure tracker capacity
    if (seenTracker.size() < pGraph.originalEdgeToPlanarEdges.size()) {
        seenTracker.resize(pGraph.originalEdgeToPlanarEdges.size(), false);
    }

    for (int pEdgeId : planarEdgesBuffer) {
        if (pGraph.hasEdge(pEdgeId)) {
            int origId = pGraph.getEdge(pEdgeId).original_edge_id;
            // O(1) Check and Insert
            if (origId >= 0 && !seenTracker[origId]) {
                seenTracker[origId] = true;
                originalEdgesBuffer.push_back(origId);
            }
        }
    }

    // O(N) Cleanup: Only reset the elements we actually touched
    for (int origId : originalEdgesBuffer) {
        seenTracker[origId] = false;
    }
}

} // namespace

std::pair<int, int> RelocationManager::normalizeEdgePair(int edgeA, int edgeB) const {
    if (edgeA <= edgeB) return {edgeA, edgeB};
    return {edgeB, edgeA};
}

// Update precomputeOriginalAdjacency to handle pre-allocation
void RelocationManager::precomputeOriginalAdjacency() {
    int maxNodeId = 0;
    m_pGraph.forEachNode([&maxNodeId](int id, const auto&) {
        if (id > maxNodeId) maxNodeId = id;
    });
    
    m_originalAdjacencyCache.assign(maxNodeId + 1, {});
    // Ensure the endpoint cache is ready for O(1) access
    m_originalEdgeEndpointsCache.assign(m_pGraph.originalEdgeToPlanarEdges.size(), {false, -1, -1});

    for (size_t origEdgeId = 0; origEdgeId < m_pGraph.originalEdgeToPlanarEdges.size(); ++origEdgeId) {
        auto ep = findOriginalEdgeEndpoints(origEdgeId);
        if (ep) {
            int u = ep->first;
            int v = ep->second;
            if (u <= maxNodeId) m_originalAdjacencyCache[u].push_back({v, origEdgeId});
            if (v <= maxNodeId) m_originalAdjacencyCache[v].push_back({u, origEdgeId});
        }
    }
}



int RelocationManager::computeActiveSideForBoundary(int currentFaceId,
                                                    const DualGraphEdge& boundary,
                                                    const DualGraph& dualGraph) const {
    if (currentFaceId < 0 || currentFaceId >= static_cast<int>(dualGraph.faces.size())) {
        return 1;
    }

    const auto& currentFacePolygon = dualGraph.faces[currentFaceId].vertices;
    if (currentFacePolygon.empty()) {
        return 1;
    }

    const double dx = boundary.x2 - boundary.x1;
    const double dy = boundary.y2 - boundary.y1;
    const double len = std::sqrt(dx * dx + dy * dy);
    if (len < kIntersectionEpsilon) return 1;

    const double midX = (boundary.x1 + boundary.x2) * 0.5;
    const double midY = (boundary.y1 + boundary.y2) * 0.5;
    const double perpX = -dy / len;
    const double perpY = dx / len;
    const double epsilon = kIntersectionEpsilon;

    const std::pair<double, double> candidateA{midX + epsilon * perpX, midY + epsilon * perpY};
    const std::pair<double, double> candidateB{midX - epsilon * perpX, midY - epsilon * perpY};

    std::pair<double, double> activePoint = candidateA;
    if (!relocation_detail::pointInPolygon(candidateA.first, candidateA.second, currentFacePolygon) &&
        relocation_detail::pointInPolygon(candidateB.first, candidateB.second, currentFacePolygon)) {
        activePoint = candidateB;
    }

    int activeSide = whichSideOfLine(
        activePoint.first,
        activePoint.second,
        boundary.x1,
        boundary.y1,
        boundary.x2,
        boundary.y2
    );

    if (activeSide == 0) return 1;
    return activeSide;
}

std::vector<int> RelocationManager::collectPathFromSource(int sourceFaceId,
                                                          int targetFaceId,
                                                          const std::vector<int>& faceParent) const {
    std::vector<int> reversedPath;

    if (sourceFaceId < 0 || targetFaceId < 0 ||
        sourceFaceId >= static_cast<int>(faceParent.size()) ||
        targetFaceId >= static_cast<int>(faceParent.size())) {
        return {};
    }

    int current = targetFaceId;
    while (current != -1) {
        reversedPath.push_back(current);
        if (current == sourceFaceId) break;
        current = faceParent[current];
    }

    if (reversedPath.empty() || reversedPath.back() != sourceFaceId) {
        return {};
    }

    std::reverse(reversedPath.begin(), reversedPath.end());
    return reversedPath;
}

std::vector<DualGraphEdge> RelocationManager::collectBoundaryPathEdges(const LocalRegionAnalysis& analysis,
                                                                       int targetFaceId) const {
    std::vector<DualGraphEdge> pathEdges;
    const auto pathFaces = collectPathFromSource(analysis.sourceFaceId, targetFaceId, analysis.faceParent);
    if (pathFaces.size() < 2) return pathEdges;

    auto findTreeEdge = [&](int a, int b) -> std::optional<DualGraphEdge> {
        for (const auto& e : analysis.dualTreeEdges) {
            if ((e.faceA == a && e.faceB == b) || (e.faceA == b && e.faceB == a)) {
                return e;
            }
        }
        return std::nullopt;
    };

    for (size_t i = 0; i + 1 < pathFaces.size(); ++i) {
        const int a = pathFaces[i];
        const int b = pathFaces[i + 1];

        auto edge = findTreeEdge(a, b);
        if (edge.has_value()) {
            pathEdges.push_back(*edge);
            continue;
        }

        for (const auto& e : analysis.dualGraph.adjacency) {
            if ((e.faceA == a && e.faceB == b) || (e.faceA == b && e.faceB == a)) {
                pathEdges.push_back(e);
                break;
            }
        }
    }

    return pathEdges;
}

// Pure O(1) Adjacency Lookup
const std::vector<std::pair<int, size_t>>& RelocationManager::collectIncidentEdgesForNode(int nodeId) const {
    // We need a static empty vector to return if the nodeId is invalid
    static const std::vector<std::pair<int, size_t>> empty_adj;

    if (nodeId < 0 || nodeId >= static_cast<int>(m_originalAdjacencyCache.size())) {
        return empty_adj;
    }
    return m_originalAdjacencyCache[nodeId];
}

// Pure math: No node lookups, no graph dependency
bool RelocationManager::intersectSegments(double x1, double y1, double x2, double y2,
                                          double x3, double y3, double x4, double y4,
                                          double& outX, double& outY) const {
    const double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (std::abs(den) < kIntersectionEpsilon) return false;

    const double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
    const double u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den;

    if (t <= kIntersectionEpsilon || t >= 1.0 - kIntersectionEpsilon || 
        u <= kIntersectionEpsilon || u >= 1.0 - kIntersectionEpsilon) {
        return false;
    }

    outX = x1 + t * (x2 - x1);
    outY = y1 + t * (y2 - y1);
    return true;
}

const std::vector<std::pair<int, size_t>>& RelocationManager::collectSelectedNeighborEdges(int variableNodeId) const {
    return collectIncidentEdgesForNode(variableNodeId);
}

void RelocationManager::collectTransitionCrossingPairs(int variableNodeId,
                                                       int currentFaceId,
                                                       const DualGraphEdge& boundary,
                                                       std::vector<std::pair<int, int>>& disappear,
                                                       std::vector<std::pair<int, int>>& appear,
                                                       const DualGraph& dualGraph) const {
    const int activeSide = computeActiveSideForBoundary(currentFaceId, boundary, dualGraph);
    const auto selectedEdges = collectSelectedNeighborEdges(variableNodeId);

    if (boundary.boundaryType == LocalSegmentType::ORIGINAL_SUBSEGMENT) {
        const int boundaryEdgeId = boundary.originalEdgeId;

        if (boundaryEdgeId < 0) {
            return;
        }

        for (const auto& [neighborId, selectedEdgeId] : selectedEdges) {
            if (selectedEdgeId < 0 || selectedEdgeId == boundaryEdgeId) {
                continue;
            }

            if (!m_pGraph.hasNode(neighborId)) {
                continue;
            }
            
            const auto& neighborNode = m_pGraph.getNode(neighborId);
            const int side = whichSideOfLine(
                neighborNode.x,
                neighborNode.y,
                boundary.x1,
                boundary.y1,
                boundary.x2,
                boundary.y2
            );

            if (side == 0) {
                continue;
            }

            const auto pair = normalizeEdgePair(selectedEdgeId, boundaryEdgeId);
            if (side == activeSide) {
                appear.push_back(pair);
            } else {
                disappear.push_back(pair);
            }
        }
        return;
    }

    if (boundary.boundaryType == LocalSegmentType::RAY) {
        const int raySourceId = boundary.raySourceNodeId;
        const int rayEmittingId = boundary.rayEmittingNodeId;

        if (raySourceId < 0 || rayEmittingId < 0) {
            return;
        }

        // Find the specific edge ID connecting the variable node to the emitting node
        int pivotEdgeId = -1;
        for (const auto& [neighborId, edgeId] : selectedEdges) {
            if (neighborId == rayEmittingId) {
                pivotEdgeId = edgeId;
                break;
            }
        }

        if (pivotEdgeId < 0) {
            return;
        }

        const auto raySourceIncident = collectIncidentEdgesForNode(raySourceId);

        for (const auto& [rayNeighborId, rayEdgeId] : raySourceIncident) {
            // Standard validity checks
            if (rayEdgeId < 0 || rayNeighborId == variableNodeId || rayEdgeId == pivotEdgeId) {
                continue;
            }

            if (!m_pGraph.hasNode(rayNeighborId)) {
                continue;
            }

            const auto& rayNeighborNode = m_pGraph.getNode(rayNeighborId);
            const int rayNeighborSide = whichSideOfLine(
                rayNeighborNode.x, rayNeighborNode.y,
                boundary.x1, boundary.y1,
                boundary.x2, boundary.y2
            );

            // Skip if the ray neighbor is collinear with the ray itself
            if (rayNeighborSide == 0) {
                continue;
            }

            const auto pair = normalizeEdgePair(pivotEdgeId, rayEdgeId);

            // LOGIC: The pivot edge sweeps from 'activeSide' to the opposite side.
            bool isAppear = (rayNeighborSide != activeSide);

            if (isAppear) {
                appear.push_back(pair);
            } else {
                disappear.push_back(pair);
            }
        }
    }
}

// Update this function in RelocationManagerCrossings.cpp
std::optional<std::pair<int, int>> RelocationManager::findOriginalEdgeEndpoints(int originalEdgeId) const {
    if (originalEdgeId < 0) return std::nullopt;
    
    // OPTIMIZATION: Pure O(1) lookup. Pre-allocation happens in precomputeOriginalAdjacency.
    if (originalEdgeId < static_cast<int>(m_originalEdgeEndpointsCache.size())) {
        const auto& cached = m_originalEdgeEndpointsCache[originalEdgeId];
        if (cached.valid) {
            return std::make_pair(cached.endpointA, cached.endpointB);
        }
    } else {
        // Fallback for safety if the graph grew unexpectedly
        m_originalEdgeEndpointsCache.resize(originalEdgeId + 100); 
    }

    if (originalEdgeId >= static_cast<int>(m_pGraph.originalEdgeToPlanarEdges.size()) ||
        m_pGraph.originalEdgeToPlanarEdges[originalEdgeId].empty()) {
        return std::nullopt;
    }

    auto traverseToOriginal = [&](int start, int previous) -> std::optional<int> {
        int current = start;
        int prev = previous;

        while (true) {
            if (!m_pGraph.hasNode(current)) return std::nullopt;
            const auto& node = m_pGraph.getNode(current);

            if (node.type == PlanarizedGraph::NodeType::ORIGINAL) {
                return current;
            }

            int next = -1;
            if (node.original_edge_1 == originalEdgeId) {
                if (node.e1_neighbor_prev == prev) next = node.e1_neighbor_next;
                else if (node.e1_neighbor_next == prev) next = node.e1_neighbor_prev;
            }
            if (next == -1 && node.original_edge_2 == originalEdgeId) {
                if (node.e2_neighbor_prev == prev) next = node.e2_neighbor_next;
                else if (node.e2_neighbor_next == prev) next = node.e2_neighbor_prev;
            }

            if (next == -1) return std::nullopt;
            prev = current;
            current = next;
        }
        return std::nullopt; // FIXED: MSVC C4715 requirement
    };

    const int anyPlanarEdgeId = m_pGraph.originalEdgeToPlanarEdges[originalEdgeId].front();
    if (!m_pGraph.hasEdge(anyPlanarEdgeId)) return std::nullopt;
    
    const auto& e = m_pGraph.getEdge(anyPlanarEdgeId);
    auto endA = traverseToOriginal(e.u_id, e.v_id);
    auto endB = traverseToOriginal(e.v_id, e.u_id);

    if (!endA.has_value() || !endB.has_value()) return std::nullopt;

    // Cache the result
    auto& cached = m_originalEdgeEndpointsCache[originalEdgeId];
    cached.valid = true;
    cached.endpointA = *endA;
    cached.endpointB = *endB;

    return std::make_pair(*endA, *endB);
}

// Refactored wrapper for existing logic that still relies on edge IDs
bool RelocationManager::intersectOriginalEdgesForMove(int edgeA, int edgeB, int variableNodeId, 
                                                      double movedX, double movedY, 
                                                      double& outX, double& outY) const {
    auto epA = findOriginalEdgeEndpoints(edgeA);
    auto epB = findOriginalEdgeEndpoints(edgeB);
    if (!epA || !epB) return false;

    auto getPoint = [&](int nodeId, double& x, double& y) {
        if (nodeId == variableNodeId) { x = movedX; y = movedY; return; }
        const auto& node = m_pGraph.getNode(nodeId); x = node.x; y = node.y;
    };

    double x1, y1, x2, y2, x3, y3, x4, y4;
    getPoint(epA->first, x1, y1); getPoint(epA->second, x2, y2);
    getPoint(epB->first, x3, y3); getPoint(epB->second, x4, y4);

    // Fast AABB check here
    double min1X = std::min(x1, x2), max1X = std::max(x1, x2);
    double min2X = std::min(x3, x4), max2X = std::max(x3, x4);
    if (max1X < min2X || min1X > max2X) return false;

    double min1Y = std::min(y1, y2), max1Y = std::max(y1, y2);
    double min2Y = std::min(y3, y4), max2Y = std::max(y3, y4);
    if (max1Y < min2Y || min1Y > max2Y) return false;

    return intersectSegments(x1, y1, x2, y2, x3, y3, x4, y4, outX, outY);
}

// 1. Refactored Counter: Zero Graph Lookups, Pure Math
int RelocationManager::countGlobalCrossingsForVariableAtPosition(
    int variableNodeId, 
    double movedX, double movedY,
    const std::unordered_set<int>& variableIncidentEdges,
    const std::vector<PrecomputedStaticEdge>& staticEdges) const {

    if (variableIncidentEdges.empty() || staticEdges.empty()) return 0;

    int count = 0;
    double outX, outY;

    for (int movedEdgeId : variableIncidentEdges) {
        auto movedEndpoints = findOriginalEdgeEndpoints(movedEdgeId);
        if (!movedEndpoints) continue;

        // Resolve moving coordinates ONCE per moved edge
        double movedAx, movedAy, movedBx, movedBy;
        if (movedEndpoints->first == variableNodeId) { movedAx = movedX; movedAy = movedY; }
        else { const auto& n = m_pGraph.getNode(movedEndpoints->first); movedAx = n.x; movedAy = n.y; }
        
        if (movedEndpoints->second == variableNodeId) { movedBx = movedX; movedBy = movedY; }
        else { const auto& n = m_pGraph.getNode(movedEndpoints->second); movedBx = n.x; movedBy = n.y; }

        // Precompute AABB of the moving edge ONCE
        const double minMx = std::min(movedAx, movedBx);
        const double maxMx = std::max(movedAx, movedBx);
        const double minMy = std::min(movedAy, movedBy);
        const double maxMy = std::max(movedAy, movedBy);

        // Inner Loop: Extremely fast contiguous array iteration
        for (const auto& staticEdge : staticEdges) {
            // Endpoint check (O(1) memory)
            if (movedEndpoints->first == staticEdge.epA || movedEndpoints->first == staticEdge.epB ||
                movedEndpoints->second == staticEdge.epA || movedEndpoints->second == staticEdge.epB) {
                continue;
            }

            // Fast AABB Check against precomputed boundaries
            if (maxMx < staticEdge.minX || minMx > staticEdge.maxX) continue;
            if (maxMy < staticEdge.minY || minMy > staticEdge.maxY) continue;

            // Pure math execution
            if (intersectSegments(movedAx, movedAy, movedBx, movedBy, 
                                  staticEdge.ax, staticEdge.ay, staticEdge.bx, staticEdge.by, 
                                  outX, outY)) {
                ++count;
            }
        }
    }
    return count;
}

int RelocationManager::evaluateGlobalCrossingDeltaForMove(
    int variableNodeId,
    double movedX,
    double movedY,
    const std::unordered_set<int>& variableIncidentEdges,
    const std::vector<PrecomputedStaticEdge>& staticEdges) const {
    
    if (!m_pGraph.hasNode(variableNodeId)) return 0;
    const auto& variableNode = m_pGraph.getNode(variableNodeId);

    // Call with the new 5th argument
    const int before = countGlobalCrossingsForVariableAtPosition(
        variableNodeId,
        variableNode.x,
        variableNode.y,
        variableIncidentEdges,
        staticEdges
    );

    // Call with the new 5th argument
    const int after = countGlobalCrossingsForVariableAtPosition(
        variableNodeId,
        movedX,
        movedY,
        variableIncidentEdges,
        staticEdges
    );

    return after - before;
}

CrossingUpdatePlan RelocationManager::buildCrossingUpdatePlan(int variableNodeId,
                                                              const LocalRegionAnalysis& analysis,
                                                              int targetFaceId) const {
    CrossingUpdatePlan plan;

    if (analysis.sourceFaceId < 0 || targetFaceId < 0 ||
        targetFaceId >= static_cast<int>(analysis.faceParent.size())) {
        return plan;
    }

    const auto boundaryPath = collectBoundaryPathEdges(analysis, targetFaceId);
    if (boundaryPath.empty()) return plan;

    std::set<std::pair<int, int>> disappearSet;
    std::set<std::pair<int, int>> appearSet;

    int currentFace = analysis.sourceFaceId;
    for (const auto& boundary : boundaryPath) {
        std::vector<std::pair<int, int>> localDisappear;
        std::vector<std::pair<int, int>> localAppear;

        collectTransitionCrossingPairs(variableNodeId,
                                       currentFace,
                                       boundary,
                                       localDisappear,
                                       localAppear,
                                       analysis.dualGraph);

        for (const auto& p : localDisappear) {
            disappearSet.insert(normalizeEdgePair(p.first, p.second));
        }
        for (const auto& p : localAppear) {
            appearSet.insert(normalizeEdgePair(p.first, p.second));
        }

        if (boundary.faceA == currentFace) currentFace = boundary.faceB;
        else if (boundary.faceB == currentFace) currentFace = boundary.faceA;
    }

    for (const auto& p : disappearSet) {
        if (appearSet.find(p) == appearSet.end()) {
            plan.disappearingPairs.push_back(p);
        }
    }

    for (const auto& p : appearSet) {
        if (disappearSet.find(p) == disappearSet.end()) {
            plan.appearingPairs.push_back(p);
        }
    }

    for (const auto& p : plan.disappearingPairs) {
        const int crossingId = m_pGraph.getCrossingNodeForPair(p.first, p.second);
        if (crossingId != -1) {
            plan.crossingNodeIdsToDestroy.push_back(crossingId);
        }
    }

    return plan;
}

// These functions are kept for header compatibility but bypassed by the Wipe and Rebuild strategy.
void RelocationManager::applyCrossingRemovals(const CrossingUpdatePlan& plan) {
    for (int crossingNodeId : plan.crossingNodeIdsToDestroy) {
        if (!m_pGraph.hasNode(crossingNodeId)) continue;
        const auto& node = m_pGraph.getNode(crossingNodeId);
        if (node.type != PlanarizedGraph::NodeType::CROSSING) continue;
        m_pGraph.destroyCrossing(crossingNodeId);
    }
}

void RelocationManager::applyCrossingInsertions(int variableNodeId,
                                                double movedX,
                                                double movedY,
                                                const CrossingUpdatePlan& plan) {
    for (const auto& p : plan.appearingPairs) {
        const int edgeA = p.first;
        const int edgeB = p.second;

        if (m_pGraph.getCrossingNodeForPair(edgeA, edgeB) != -1) continue;

        double intX, intY;
        if (intersectOriginalEdgesForMove(edgeA, edgeB, variableNodeId, movedX, movedY, intX, intY)) {
            m_pGraph.createCrossing(edgeA, edgeB, intX, intY);
        }
    }
}

void RelocationManager::reconcileCrossingsForMovedEdges(
    const std::unordered_set<int>& movedOriginalEdges,
    const CrossingUpdatePlan& plan,
    int variableNodeId,
    double movedX, double movedY
) {
    if (movedOriginalEdges.empty()) return;

    const auto& grid = m_pGraph.getGrid();
    if (grid.getNumCells() == 0) return;

    // 1. Verify appearing pairs from the plan
    for (const auto& planPair : plan.appearingPairs) {
        double intX, intY;
        if (intersectOriginalEdgesForMove(planPair.first, planPair.second, variableNodeId, movedX, movedY, intX, intY)) {
            if (m_pGraph.getCrossingNodeForPair(planPair.first, planPair.second) == -1) {
                m_pGraph.createCrossing(planPair.first, planPair.second, intX, intY);
            }
        }
    }

    // 2. Search for additional new/remaining intersections
    for (int movedEdgeId : movedOriginalEdges) {
        auto ep = findOriginalEdgeEndpoints(movedEdgeId);
        if (!ep) continue;

        // NEW: Hoist the moving edge coordinate resolution outside the inner loop!
        double movedAx, movedAy, movedBx, movedBy;
        if (ep->first == variableNodeId) { movedAx = movedX; movedAy = movedY; }
        else { const auto& nA = m_pGraph.getNode(ep->first); movedAx = nA.x; movedAy = nA.y; }

        if (ep->second == variableNodeId) { movedBx = movedX; movedBy = movedY; }
        else { const auto& nB = m_pGraph.getNode(ep->second); movedBx = nB.x; movedBy = nB.y; }

        // NEW: Precompute the AABB of the moving edge
        const double minMx = std::min(movedAx, movedBx);
        const double maxMx = std::max(movedAx, movedBx);
        const double minMy = std::min(movedAy, movedBy);
        const double maxMy = std::max(movedAy, movedBy);

        collectCandidateOriginalEdgesForSegment(
            m_pGraph, grid, movedAx, movedAy, movedBx, movedBy,
            m_tempPlanarEdgesBuffer, m_tempOriginalEdgesBuffer, m_seenOriginalEdgesTracker);

        for (int otherEdgeId : m_tempOriginalEdgesBuffer) {
            if (movedEdgeId == otherEdgeId) continue;
            if (movedOriginalEdges.count(otherEdgeId) > 0 && movedEdgeId > otherEdgeId) continue;

            auto otherEp = findOriginalEdgeEndpoints(otherEdgeId);
            if (!otherEp) continue;

            if (ep->first == otherEp->first || ep->first == otherEp->second ||
                ep->second == otherEp->first || ep->second == otherEp->second) continue;

            // NEW: Resolve static edge coordinates and run fast AABB inline
            const auto& nC = m_pGraph.getNode(otherEp->first);
            const auto& nD = m_pGraph.getNode(otherEp->second);

            const double minOx = std::min(nC.x, nD.x);
            const double maxOx = std::max(nC.x, nD.x);
            if (maxMx < minOx || minMx > maxOx) continue;

            const double minOy = std::min(nC.y, nD.y);
            const double maxOy = std::max(nC.y, nD.y);
            if (maxMy < minOy || minMy > maxOy) continue;

            // NEW: Bypass intersectOriginalEdgesForMove overhead entirely
            double otherIntX, otherIntY;
            if (intersectSegments(movedAx, movedAy, movedBx, movedBy, nC.x, nC.y, nD.x, nD.y, otherIntX, otherIntY)) {
                const int existingId = m_pGraph.getCrossingNodeForPair(movedEdgeId, otherEdgeId);
                if (existingId == -1) {
                    m_pGraph.createCrossing(movedEdgeId, otherEdgeId, otherIntX, otherIntY);
                }
            }
        }
    }
}

RelocationStepResult RelocationManager::performRelocationStep() {
    RelocationStepResult result;

    const int variableNodeId = selectVariableNode();
    if (variableNodeId == -1) return result;

    result.valid = true;
    result.variableNodeId = variableNodeId;

    const RegionOfInterest roi = calculateROI(variableNodeId);
    const LocalRegionAnalysis analysis = analyzeLocalRegions(roi, variableNodeId);
    result.sourceFaceId = analysis.sourceFaceId;

    std::unordered_set<int> movedOriginalEdges;
    for (const auto& [neighborId, originalEdgeId] : collectSelectedNeighborEdges(variableNodeId)) {
        (void)neighborId;
        if (originalEdgeId >= 0) movedOriginalEdges.insert(originalEdgeId);
    }
    if (movedOriginalEdges.empty()) {
        return result;
    }

    std::random_device rd;
    std::mt19937 rng(rd());

    auto chosenTarget = chooseTargetFace(analysis, rng);
    if (!chosenTarget.has_value()) {
        return result;
    }

    const int targetFaceId = chosenTarget->first;
    const double targetWeight = chosenTarget->second;

    if (targetFaceId < 0 || targetFaceId >= static_cast<int>(analysis.dualGraph.faces.size())) {
        return result;
    }

    const auto targetPoint = chooseInteriorPointInFace(analysis.dualGraph.faces[targetFaceId]);
    if (!targetPoint.has_value()) {
        return result;
    }

    const double movedX = targetPoint->first;
    const double movedY = targetPoint->second;

    CrossingUpdatePlan plan = buildCrossingUpdatePlan(variableNodeId, analysis, targetFaceId);

    // =========================================================================
    // EXECUTE GRAPH UPDATES: The "Wipe and Rebuild" Strategy
    // =========================================================================
    // 1. WIPE: Destroy all crossings on moved edges to guarantee pristine straight lines.
    std::unordered_set<int> crossingsToWipe;
    for (int movedEdgeId : movedOriginalEdges) {
        if (movedEdgeId >= 0 && movedEdgeId < static_cast<int>(m_pGraph.originalEdgeToPlanarEdges.size())) {
            for (int pEdgeId : m_pGraph.originalEdgeToPlanarEdges[movedEdgeId]) {
                if (!m_pGraph.hasEdge(pEdgeId)) continue;
                int u = m_pGraph.getEdge(pEdgeId).u_id;
                int v = m_pGraph.getEdge(pEdgeId).v_id;
                if (m_pGraph.getNode(u).type == PlanarizedGraph::NodeType::CROSSING) crossingsToWipe.insert(u);
                if (m_pGraph.getNode(v).type == PlanarizedGraph::NodeType::CROSSING) crossingsToWipe.insert(v);
            }
        }
    }
    for (int cId : crossingsToWipe) {
        m_pGraph.destroyCrossing(cId);
    }

    // 2. MOVE NODE: Lines are now perfectly straight segments going to the new coordinate.
    m_pGraph.updateNodePosition(variableNodeId, movedX, movedY);

    // 3. REBUILD: Let Reconcile safely reconstruct all necessary crossings on straight geometry.
    // (applyCrossingRemovals and applyCrossingInsertions are intentionally bypassed)
    reconcileCrossingsForMovedEdges(movedOriginalEdges, plan, variableNodeId, movedX, movedY);

    result.moved = true;
    result.usedExplorationMove = false;
    result.targetFaceId = targetFaceId;
    result.targetWeight = targetWeight;
    result.targetGlobalDelta =
        (targetFaceId >= 0 && targetFaceId < static_cast<int>(analysis.faceGlobalCrossingDelta.size()))
            ? analysis.faceGlobalCrossingDelta[targetFaceId]
            : 0;
    return result;
}