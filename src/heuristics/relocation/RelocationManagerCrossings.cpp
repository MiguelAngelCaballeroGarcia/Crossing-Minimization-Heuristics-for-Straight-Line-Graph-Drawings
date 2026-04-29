#include "../RelocationManager.hpp"
#include "RelocationHelpers.hpp"

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

namespace {

constexpr double kIntersectionEpsilon = 1e-9;

inline std::optional<std::pair<double, double>> pureMathIntersect(
    double x1, double y1, double x2, double y2,
    double x3, double y3, double x4, double y4) 
{
    // Heavy Math
    const double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (std::abs(den) < kIntersectionEpsilon) return std::nullopt;

    const double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
    const double u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den;

    const double eps = kIntersectionEpsilon;
    if (t <= eps || t >= 1.0 - eps || u <= eps || u >= 1.0 - eps) {
        return std::nullopt;
    }

    return std::make_pair(x1 + t * (x2 - x1), y1 + t * (y2 - y1));
}

// MODIFICATION: Pass output container by reference to prevent heap allocation churn
void collectCandidateEdgesForSegment(const SpatialGrid& grid,
                                     double x1, double y1,
                                     double x2, double y2,
                                     std::unordered_set<int>& outCandidateEdges) {
    if (grid.getNumCells() == 0) return;

    const double minX = grid.getMinX();
    const double minY = grid.getMinY();
    const double cellWidth = grid.getCellWidth();
    const double cellHeight = grid.getCellHeight();

    int minCol = static_cast<int>(std::floor((std::min(x1, x2) - minX) / cellWidth));
    int maxCol = static_cast<int>(std::floor((std::max(x1, x2) - minX) / cellWidth));
    int minRow = static_cast<int>(std::floor((std::min(y1, y2) - minY) / cellHeight));
    int maxRow = static_cast<int>(std::floor((std::max(y1, y2) - minY) / cellHeight));

    minCol = std::clamp(minCol, 0, grid.getNumCellsX() - 1);
    maxCol = std::clamp(maxCol, 0, grid.getNumCellsX() - 1);
    minRow = std::clamp(minRow, 0, grid.getNumCellsY() - 1);
    maxRow = std::clamp(maxRow, 0, grid.getNumCellsY() - 1);

    for (int col = minCol; col <= maxCol; ++col) {
        for (int row = minRow; row <= maxRow; ++row) {
            const int cellIndex = relocation_detail::toCellIndex(col, row, grid);
            const auto& cell = grid.getCell(cellIndex);
            outCandidateEdges.insert(cell.edgeIndices.begin(), cell.edgeIndices.end());
        }
    }
}

// Add to RelocationManagerCrossings.cpp anonymous namespace or as a private member
inline uint64_t makeOriginalEdgePairKey(int edgeA, int edgeB) {
    int a = std::min(edgeA, edgeB);
    int b = std::max(edgeA, edgeB);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b);
}

// MODIFICATION: Updated to match the pass-by-reference signature
void collectCandidateOriginalEdgesForSegment(const PlanarizedGraph& pGraph,
                                             const SpatialGrid& grid,
                                             double x1, double y1,
                                             double x2, double y2,
                                             std::unordered_set<int>& outCandidateOriginalEdges) {
    std::unordered_set<int> candidatePlanarEdges;
    collectCandidateEdgesForSegment(grid, x1, y1, x2, y2, candidatePlanarEdges);

    for (int planarEdgeId : candidatePlanarEdges) {
        if (!pGraph.hasEdge(planarEdgeId)) continue;
        outCandidateOriginalEdges.insert(pGraph.getEdge(planarEdgeId).original_edge_id);
    }
}

} // namespace

// Forward declaration for helper implemented later in this file
void printCrossingDifferenceReport(const PlanarizedGraph& pGraph);

std::pair<int, int> RelocationManager::normalizeEdgePair(int edgeA, int edgeB) const {
    if (edgeA <= edgeB) return {edgeA, edgeB};
    return {edgeB, edgeA};
}

void RelocationManager::resetStepCrossingCaches() const {
    m_pairToOriginalEdgeCache.clear();
    m_pairToOriginalEdgeCacheValid = false;
    m_originalEdgeEndpointsCache.clear();
}

void RelocationManager::ensurePairToOriginalEdgeCache() const {
    if (m_pairToOriginalEdgeCacheValid) {
        return;
    }

    m_pairToOriginalEdgeCache.clear();
    m_pairToOriginalEdgeCache.reserve(m_pGraph.edges.size());
    m_pGraph.forEachEdge([&](int, const PlanarizedGraph::PlanarEdge& edge) {
        const int a = std::min(edge.u_id, edge.v_id);
        const int b = std::max(edge.u_id, edge.v_id);
        const long long key = (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);
        m_pairToOriginalEdgeCache[key] = edge.original_edge_id;
    });

    m_pairToOriginalEdgeCacheValid = true;
}

std::vector<std::pair<int, int>> RelocationManager::collectIncidentEdgesForNode(int nodeId) const {
    std::vector<std::pair<int, int>> incident;

    if (!m_pGraph.hasNode(nodeId)) return incident;
    const auto& startNode = m_pGraph.getNode(nodeId);
    if (startNode.type != PlanarizedGraph::NodeType::ORIGINAL) return incident;

    ensurePairToOriginalEdgeCache();

    std::unordered_set<int> seenOriginalEdges;
    for (int firstHopId : startNode.adjacent_planar_nodes) {
        const int a = std::min(nodeId, firstHopId);
        const int b = std::max(nodeId, firstHopId);
        const long long key = (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);

        auto edgeIt = m_pairToOriginalEdgeCache.find(key);
        if (edgeIt == m_pairToOriginalEdgeCache.end()) continue;
        const int originalEdgeId = edgeIt->second;
        if (!seenOriginalEdges.insert(originalEdgeId).second) continue;

        int current = firstHopId;
        int previous = nodeId;

        while (true) {
            if (!m_pGraph.hasNode(current)) break;

            const auto& currNode = m_pGraph.getNode(current);
            if (currNode.type == PlanarizedGraph::NodeType::ORIGINAL) {
                incident.push_back({current, originalEdgeId});
                break;
            }

            int next = -1;
            if (currNode.original_edge_1 == originalEdgeId) {
                if (currNode.e1_neighbor_prev == previous) next = currNode.e1_neighbor_next;
                else if (currNode.e1_neighbor_next == previous) next = currNode.e1_neighbor_prev;
            }
            if (next == -1 && currNode.original_edge_2 == originalEdgeId) {
                if (currNode.e2_neighbor_prev == previous) next = currNode.e2_neighbor_next;
                else if (currNode.e2_neighbor_next == previous) next = currNode.e2_neighbor_prev;
            }

            if (next == -1) break;
            previous = current;
            current = next;
        }
    }

    return incident;
}

std::vector<std::pair<int, int>> RelocationManager::collectSelectedNeighborEdges(int variableNodeId) const {
    return collectIncidentEdgesForNode(variableNodeId);
}

std::optional<std::pair<int, int>> RelocationManager::findOriginalEdgeEndpoints(int originalEdgeId) const {
    auto cached = m_originalEdgeEndpointsCache.find(originalEdgeId);
    if (cached != m_originalEdgeEndpointsCache.end()) {
        return cached->second;
    }

    if (originalEdgeId < 0 ||
        originalEdgeId >= static_cast<int>(m_pGraph.originalEdgeToPlanarEdges.size()) ||
        m_pGraph.originalEdgeToPlanarEdges[originalEdgeId].empty()) {
        m_originalEdgeEndpointsCache[originalEdgeId] = std::nullopt;
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
    };

    const int anyPlanarEdgeId = m_pGraph.originalEdgeToPlanarEdges[originalEdgeId].front();
    if (!m_pGraph.hasEdge(anyPlanarEdgeId)) {
        m_originalEdgeEndpointsCache[originalEdgeId] = std::nullopt;
        return std::nullopt;
    }
    const auto& e = m_pGraph.getEdge(anyPlanarEdgeId);

    const int a = e.u_id;
    const int b = e.v_id;
    auto endA = traverseToOriginal(a, b);
    auto endB = traverseToOriginal(b, a);

    if (!endA.has_value() || !endB.has_value()) {
        m_originalEdgeEndpointsCache[originalEdgeId] = std::nullopt;
        return std::nullopt;
    }

    const auto endpoints = std::make_pair(*endA, *endB);
    m_originalEdgeEndpointsCache[originalEdgeId] = endpoints;
    return endpoints;
}

std::optional<std::pair<double, double>> RelocationManager::intersectOriginalEdgesForMove(int edgeA,
                                                                                          int edgeB,
                                                                                          int variableNodeId,
                                                                                          double movedX,
                                                                                          double movedY) const {
    auto epA = findOriginalEdgeEndpoints(edgeA);
    auto epB = findOriginalEdgeEndpoints(edgeB);
    if (!epA.has_value() || !epB.has_value()) return std::nullopt;

    auto getPointForEdgeEndpoint = [&](int nodeId) -> std::optional<std::pair<double, double>> {
        if (nodeId == variableNodeId) {
            return std::make_pair(movedX, movedY);
        }
        if (!m_pGraph.hasNode(nodeId)) return std::nullopt;
        const auto& node = m_pGraph.getNode(nodeId);
        return std::make_pair(node.x, node.y);
    };

    auto a1 = getPointForEdgeEndpoint(epA->first);
    auto a2 = getPointForEdgeEndpoint(epA->second);
    auto b1 = getPointForEdgeEndpoint(epB->first);
    auto b2 = getPointForEdgeEndpoint(epB->second);
    if (!a1.has_value() || !a2.has_value() || !b1.has_value() || !b2.has_value()) return std::nullopt;

    const double x1 = a1->first;
    const double y1 = a1->second;
    const double x2 = a2->first;
    const double y2 = a2->second;
    const double x3 = b1->first;
    const double y3 = b1->second;
    const double x4 = b2->first;
    const double y4 = b2->second;

    const double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (std::abs(den) < kIntersectionEpsilon) return std::nullopt;

    const double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
    const double u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den;

    const double eps = kIntersectionEpsilon;
    if (t <= eps || t >= 1.0 - eps || u <= eps || u >= 1.0 - eps) {
        return std::nullopt;
    }

    return std::make_pair(x1 + t * (x2 - x1), y1 + t * (y2 - y1));
}

int RelocationManager::countGlobalCrossingsForVariableAtPosition(
    int variableNodeId,
    double movedX,
    double movedY,
    const std::unordered_set<int>& variableIncidentEdges) const {
    
    if (variableIncidentEdges.empty()) return 0;

    const auto& grid = m_pGraph.getGrid();
    if (grid.getNumCells() == 0) return 0;

    int count = 0;
    
    std::unordered_set<uint64_t> visitedPairs;
    visitedPairs.reserve(variableIncidentEdges.size() * 10);

    // MODIFICATION: Reuse buffers passed by reference
    std::unordered_set<int> planarCandidates;
    std::unordered_set<int> candidateOriginalEdges;
    planarCandidates.reserve(256);
    candidateOriginalEdges.reserve(256);

    for (int movedEdgeId : variableIncidentEdges) {
        auto movedEndpoints = findOriginalEdgeEndpoints(movedEdgeId);
        if (!movedEndpoints.has_value()) continue;

        const int u1 = movedEndpoints->first;
        const int v1 = movedEndpoints->second;

        const double x1 = (u1 == variableNodeId) ? movedX : m_pGraph.getNode(u1).x;
        const double y1 = (u1 == variableNodeId) ? movedY : m_pGraph.getNode(u1).y;
        const double x2 = (v1 == variableNodeId) ? movedX : m_pGraph.getNode(v1).x;
        const double y2 = (v1 == variableNodeId) ? movedY : m_pGraph.getNode(v1).y;

        // MODIFICATION: Pre-calculate the bounding box for the moved segment
        const double minX1 = std::min(x1, x2);
        const double maxX1 = std::max(x1, x2);
        const double minY1 = std::min(y1, y2);
        const double maxY1 = std::max(y1, y2);

        planarCandidates.clear();
        candidateOriginalEdges.clear();
        
        collectCandidateEdgesForSegment(grid, x1, y1, x2, y2, planarCandidates);
        for (int pEdgeId : planarCandidates) {
            if (m_pGraph.hasEdge(pEdgeId)) {
                candidateOriginalEdges.insert(m_pGraph.getEdge(pEdgeId).original_edge_id);
            }
        }

        for (int otherEdgeId : candidateOriginalEdges) {
            if (otherEdgeId == movedEdgeId || variableIncidentEdges.count(otherEdgeId) > 0) {
                continue;
            }

            uint64_t pairKey = makeOriginalEdgePairKey(movedEdgeId, otherEdgeId);
            if (!visitedPairs.insert(pairKey).second) continue;

            auto otherEndpoints = findOriginalEdgeEndpoints(otherEdgeId);
            if (!otherEndpoints.has_value()) continue;

            const int u2 = otherEndpoints->first;
            const int v2 = otherEndpoints->second;
            if (u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2) {
                continue;
            }

            const double x3 = m_pGraph.getNode(u2).x;
            const double y3 = m_pGraph.getNode(u2).y;
            const double x4 = m_pGraph.getNode(v2).x;
            const double y4 = m_pGraph.getNode(v2).y;

            // MODIFICATION: Use pre-calculated AABB bounds to speed up early rejection
            if (maxX1 < std::min(x3, x4) || minX1 > std::max(x3, x4) ||
                maxY1 < std::min(y3, y4) || minY1 > std::max(y3, y4)) {
                continue;
            }

            const double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            if (std::abs(den) < kIntersectionEpsilon) continue;

            const double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
            const double u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den;

            if (t > kIntersectionEpsilon && t < 1.0 - kIntersectionEpsilon && 
                u > kIntersectionEpsilon && u < 1.0 - kIntersectionEpsilon) {
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
    const std::unordered_set<int>& variableIncidentEdges) const {
    if (!m_pGraph.hasNode(variableNodeId)) return 0;
    const auto& variableNode = m_pGraph.getNode(variableNodeId);

    const double currentX = variableNode.x;
    const double currentY = variableNode.y;

    const int before = countGlobalCrossingsForVariableAtPosition(
        variableNodeId,
        currentX,
        currentY,
        variableIncidentEdges
    );

    const int after = countGlobalCrossingsForVariableAtPosition(
        variableNodeId,
        movedX,
        movedY,
        variableIncidentEdges
    );

    return after - before;
}

RelocationStepResult RelocationManager::performRelocationStep() {
    RelocationStepResult result;

    const int variableNodeId = selectVariableNode();
    if (variableNodeId == -1) return result;

    result.valid = true;
    result.variableNodeId = variableNodeId;

    const RegionOfInterest roi = calculateROI(variableNodeId);
    
    // Phase 1: Dual-Graph Evaluation
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
    if (targetFaceId < 0 || targetFaceId >= static_cast<int>(analysis.dualGraph.faces.size())) return result;

    const auto targetPoint = chooseInteriorPointInFace(analysis.dualGraph.faces[targetFaceId]);
    if (!targetPoint.has_value()) return result;

    const double movedX = targetPoint->first;
    const double movedY = targetPoint->second;

    std::unordered_set<int> crossingsToDestroy;
    for (int movedEdgeId : movedOriginalEdges) {
        if (movedEdgeId < 0 || movedEdgeId >= static_cast<int>(m_pGraph.originalEdgeToPlanarEdges.size())) continue;
        
        for (int pEdgeId : m_pGraph.originalEdgeToPlanarEdges[movedEdgeId]) {
            if (!m_pGraph.hasEdge(pEdgeId)) continue;
            const auto& pe = m_pGraph.getEdge(pEdgeId);
            
            auto checkNode = [&](int nodeId) {
                if (m_pGraph.hasNode(nodeId) && m_pGraph.getNode(nodeId).type == PlanarizedGraph::NodeType::CROSSING) {
                    const auto& cNode = m_pGraph.getNode(nodeId);
                    if (cNode.original_edge_1 == movedEdgeId || cNode.original_edge_2 == movedEdgeId) {
                        crossingsToDestroy.insert(nodeId);
                    }
                }
            };
            checkNode(pe.u_id);
            checkNode(pe.v_id);
        }
    }
    
    if (!crossingsToDestroy.empty()) {
        for (int crossingId : crossingsToDestroy) {
            m_pGraph.destroyCrossing(crossingId);
        }
    }

    m_pGraph.updateNodePosition(variableNodeId, movedX, movedY);

    int createdCount = 0;
    
    std::unordered_set<std::uint64_t> processedPairs;
    processedPairs.reserve(movedOriginalEdges.size() * 10); 
    
    const auto& grid = m_pGraph.getGrid();
    
    // MODIFICATION: Reusable sets passed by reference
    std::unordered_set<int> planarCandidates;
    std::unordered_set<int> candidateOtherEdges;
    planarCandidates.reserve(256);
    candidateOtherEdges.reserve(256);

    for (int movedEdgeId : movedOriginalEdges) {
        auto movedEndpoints = findOriginalEdgeEndpoints(movedEdgeId);
        if (!movedEndpoints.has_value()) continue;
        
        const double aX = m_pGraph.getNode(movedEndpoints->first).x;
        const double aY = m_pGraph.getNode(movedEndpoints->first).y;
        const double bX = m_pGraph.getNode(movedEndpoints->second).x;
        const double bY = m_pGraph.getNode(movedEndpoints->second).y;
        
        // MODIFICATION: Pre-calculate bounds
        const double minA = std::min(aX, bX);
        const double maxA = std::max(aX, bX);
        const double minB = std::min(aY, bY);
        const double maxB = std::max(aY, bY);

        planarCandidates.clear();
        candidateOtherEdges.clear();
        collectCandidateEdgesForSegment(grid, aX, aY, bX, bY, planarCandidates);
        
        for (int pEdgeId : planarCandidates) {
            if (m_pGraph.hasEdge(pEdgeId)) {
                candidateOtherEdges.insert(m_pGraph.getEdge(pEdgeId).original_edge_id);
            }
        }
        
        for (int otherEdgeId : candidateOtherEdges) {
            if (movedEdgeId == otherEdgeId) continue;
            if (movedOriginalEdges.count(otherEdgeId) > 0 && movedEdgeId > otherEdgeId) continue; 
            
            std::uint64_t key = makeOriginalEdgePairKey(movedEdgeId, otherEdgeId);
            if (!processedPairs.insert(key).second) continue;
            
            auto otherEndpoints = findOriginalEdgeEndpoints(otherEdgeId);
            if (!otherEndpoints.has_value()) continue;
            
            if (movedEndpoints->first == otherEndpoints->first || movedEndpoints->first == otherEndpoints->second ||
                movedEndpoints->second == otherEndpoints->first || movedEndpoints->second == otherEndpoints->second) {
                continue;
            }

            const double cX = m_pGraph.getNode(otherEndpoints->first).x;
            const double cY = m_pGraph.getNode(otherEndpoints->first).y;
            const double dX = m_pGraph.getNode(otherEndpoints->second).x;
            const double dY = m_pGraph.getNode(otherEndpoints->second).y;
            
            // MODIFICATION: Apply AABB filter to fast-fail intersection checks
            if (maxA < std::min(cX, dX) || minA > std::max(cX, dX) || 
                maxB < std::min(cY, dY) || minB > std::max(cY, dY)) {
                continue;
            }

            auto intersection = pureMathIntersect(aX, aY, bX, bY, cX, cY, dX, dY);
            
            if (intersection.has_value()) {
                m_pGraph.createCrossing(movedEdgeId, otherEdgeId, intersection->first, intersection->second);
                createdCount++;
            }
        }
    }

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