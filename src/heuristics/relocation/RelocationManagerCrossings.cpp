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

namespace {

constexpr double kIntersectionEpsilon = 1e-9;

std::unordered_set<int> collectCandidateEdgesForSegment(const SpatialGrid& grid,
                                                        double x1,
                                                        double y1,
                                                        double x2,
                                                        double y2) {
    std::unordered_set<int> candidateEdges;
    if (grid.getNumCells() == 0) return candidateEdges;

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
            candidateEdges.insert(cell.edgeIndices.begin(), cell.edgeIndices.end());
        }
    }

    return candidateEdges;
}

} // namespace

std::pair<int, int> RelocationManager::normalizeEdgePair(int edgeA, int edgeB) const {
    if (edgeA <= edgeB) return {edgeA, edgeB};
    return {edgeB, edgeA};
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

std::vector<std::pair<int, int>> RelocationManager::collectIncidentEdgesForNode(int nodeId) const {
    std::vector<std::pair<int, int>> incident;

    if (!m_pGraph.hasNode(nodeId)) return incident;
    const auto& startNode = m_pGraph.getNode(nodeId);
    if (startNode.type != PlanarizedGraph::NodeType::ORIGINAL) return incident;

    std::unordered_map<long long, int> pairToOriginalEdge;
    m_pGraph.forEachEdge([&](int, const PlanarizedGraph::PlanarEdge& edge) {
        const int a = std::min(edge.u_id, edge.v_id);
        const int b = std::max(edge.u_id, edge.v_id);
        const long long key = (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);
        pairToOriginalEdge[key] = edge.original_edge_id;
    });

    std::unordered_set<int> seenOriginalEdges;
    for (int firstHopId : startNode.adjacent_planar_nodes) {
        const int a = std::min(nodeId, firstHopId);
        const int b = std::max(nodeId, firstHopId);
        const long long key = (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);

        auto edgeIt = pairToOriginalEdge.find(key);
        if (edgeIt == pairToOriginalEdge.end()) continue;
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
        if (boundaryEdgeId < 0) return;

        for (const auto& [neighborId, selectedEdgeId] : selectedEdges) {
            if (selectedEdgeId < 0 || selectedEdgeId == boundaryEdgeId) continue;

            if (!m_pGraph.hasNode(neighborId)) continue;
            const auto& neighborNode = m_pGraph.getNode(neighborId);

            const int side = whichSideOfLine(
                neighborNode.x,
                neighborNode.y,
                boundary.x1,
                boundary.y1,
                boundary.x2,
                boundary.y2
            );
            if (side == 0) continue;

            const auto pair = normalizeEdgePair(selectedEdgeId, boundaryEdgeId);
            if (side == activeSide) appear.push_back(pair);
            else disappear.push_back(pair);
        }
        return;
    }

    if (boundary.boundaryType == LocalSegmentType::RAY) {
        const int raySourceId = boundary.raySourceNodeId;
        if (raySourceId < 0) return;

        const auto raySourceIncident = collectIncidentEdgesForNode(raySourceId);
        for (const auto& [selectedNeighborId, selectedEdgeId] : selectedEdges) {
            if (selectedEdgeId < 0) continue;

            if (!m_pGraph.hasNode(selectedNeighborId)) continue;
            const auto& selectedNode = m_pGraph.getNode(selectedNeighborId);
            const int selectedSide = whichSideOfLine(
                selectedNode.x,
                selectedNode.y,
                boundary.x1,
                boundary.y1,
                boundary.x2,
                boundary.y2
            );
            if (selectedSide == 0) continue;

            for (const auto& [rayNeighborId, rayEdgeId] : raySourceIncident) {
                if (rayEdgeId < 0 || rayNeighborId == variableNodeId) continue;
                if (rayEdgeId == selectedEdgeId) continue;

                if (!m_pGraph.hasNode(rayNeighborId)) continue;
                const auto& rayNeighborNode = m_pGraph.getNode(rayNeighborId);

                const int rayNeighborSide = whichSideOfLine(
                    rayNeighborNode.x,
                    rayNeighborNode.y,
                    boundary.x1,
                    boundary.y1,
                    boundary.x2,
                    boundary.y2
                );
                if (rayNeighborSide == 0) continue;

                const auto pair = normalizeEdgePair(selectedEdgeId, rayEdgeId);
                if (selectedSide == activeSide) appear.push_back(pair);
                else disappear.push_back(pair);
            }
        }
    }
}

std::optional<std::pair<int, int>> RelocationManager::findOriginalEdgeEndpoints(int originalEdgeId) const {
    if (originalEdgeId < 0 ||
        originalEdgeId >= static_cast<int>(m_pGraph.originalEdgeToPlanarEdges.size()) ||
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
    };

    const int anyPlanarEdgeId = m_pGraph.originalEdgeToPlanarEdges[originalEdgeId].front();
    if (!m_pGraph.hasEdge(anyPlanarEdgeId)) return std::nullopt;
    const auto& e = m_pGraph.getEdge(anyPlanarEdgeId);

    const int a = e.u_id;
    const int b = e.v_id;
    auto endA = traverseToOriginal(a, b);
    auto endB = traverseToOriginal(b, a);

    if (!endA.has_value() || !endB.has_value()) return std::nullopt;
    return std::make_pair(*endA, *endB);
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

int RelocationManager::findExistingCrossingNodeForPair(int edgeA, int edgeB) const {
    const auto target = normalizeEdgePair(edgeA, edgeB);

    int found = -1;
    m_pGraph.forEachNode([&](int nodeId, const PlanarizedGraph::PlanarNode& node) {
        if (found != -1) return;
        if (node.type != PlanarizedGraph::NodeType::CROSSING) return;

        const auto crossingPair = normalizeEdgePair(node.original_edge_1, node.original_edge_2);
        if (crossingPair == target) {
            found = nodeId;
        }
    });

    return found;
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
    std::set<std::pair<int, int>> visitedPairs;

    for (int movedEdgeId : variableIncidentEdges) {
        auto movedEndpoints = findOriginalEdgeEndpoints(movedEdgeId);
        if (!movedEndpoints.has_value()) continue;

        auto movedPointForNode = [&](int nodeId) -> std::optional<std::pair<double, double>> {
            if (nodeId == variableNodeId) {
                return std::make_pair(movedX, movedY);
            }
            if (!m_pGraph.hasNode(nodeId)) return std::nullopt;
            const auto& node = m_pGraph.getNode(nodeId);
            return std::make_pair(node.x, node.y);
        };

        auto movedA = movedPointForNode(movedEndpoints->first);
        auto movedB = movedPointForNode(movedEndpoints->second);
        if (!movedA.has_value() || !movedB.has_value()) continue;

        const auto candidateEdges = collectCandidateEdgesForSegment(
            grid,
            movedA->first,
            movedA->second,
            movedB->first,
            movedB->second
        );

        for (int otherEdgeId : candidateEdges) {
            if (otherEdgeId == movedEdgeId) continue;
            if (variableIncidentEdges.count(otherEdgeId) > 0) continue;

            const auto pair = normalizeEdgePair(movedEdgeId, otherEdgeId);
            if (!visitedPairs.insert(pair).second) continue;

            auto otherEndpoints = findOriginalEdgeEndpoints(otherEdgeId);
            if (!otherEndpoints.has_value()) continue;

            const bool sharesEndpoint =
                movedEndpoints->first == otherEndpoints->first ||
                movedEndpoints->first == otherEndpoints->second ||
                movedEndpoints->second == otherEndpoints->first ||
                movedEndpoints->second == otherEndpoints->second;
            if (sharesEndpoint) continue;

            auto intersection = intersectOriginalEdgesForMove(
                movedEdgeId,
                otherEdgeId,
                variableNodeId,
                movedX,
                movedY
            );

            if (intersection.has_value()) {
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

        for (const auto& p : localDisappear) disappearSet.insert(normalizeEdgePair(p.first, p.second));
        for (const auto& p : localAppear) appearSet.insert(normalizeEdgePair(p.first, p.second));

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
        const int crossingId = findExistingCrossingNodeForPair(p.first, p.second);
        if (crossingId != -1) {
            plan.crossingNodeIdsToDestroy.push_back(crossingId);
        }
    }

    return plan;
}

void RelocationManager::applyCrossingRemovals(const CrossingUpdatePlan& plan) {
    for (int crossingNodeId : plan.crossingNodeIdsToDestroy) {
        if (!m_pGraph.hasNode(crossingNodeId)) continue;
        if (m_pGraph.getNode(crossingNodeId).type != PlanarizedGraph::NodeType::CROSSING) continue;
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

        if (findExistingCrossingNodeForPair(edgeA, edgeB) != -1) {
            continue;
        }

        const auto intersection = intersectOriginalEdgesForMove(edgeA, edgeB, variableNodeId, movedX, movedY);
        if (!intersection.has_value()) continue;

        m_pGraph.createCrossing(edgeA, edgeB, intersection->first, intersection->second);
    }
}

void RelocationManager::reconcileCrossingsForMovedEdges(const std::unordered_set<int>& movedOriginalEdges) {
    if (movedOriginalEdges.empty()) return;

    std::vector<int> allOriginalEdges;
    allOriginalEdges.reserve(m_pGraph.originalEdgeToPlanarEdges.size());
    for (int origEdgeId = 0; origEdgeId < static_cast<int>(m_pGraph.originalEdgeToPlanarEdges.size()); ++origEdgeId) {
        if (m_pGraph.originalEdgeToPlanarEdges[origEdgeId].empty()) continue;
        allOriginalEdges.push_back(origEdgeId);
    }

    std::set<std::pair<int, int>> processed;

    for (int movedEdgeId : movedOriginalEdges) {
        for (int otherEdgeId : allOriginalEdges) {
            if (movedEdgeId == otherEdgeId) continue;

            const auto pair = normalizeEdgePair(movedEdgeId, otherEdgeId);
            if (!processed.insert(pair).second) continue;

            auto epA = findOriginalEdgeEndpoints(pair.first);
            auto epB = findOriginalEdgeEndpoints(pair.second);
            if (!epA.has_value() || !epB.has_value()) continue;

            const bool sharesEndpoint =
                epA->first == epB->first || epA->first == epB->second ||
                epA->second == epB->first || epA->second == epB->second;
            if (sharesEndpoint) {
                const int existingId = findExistingCrossingNodeForPair(pair.first, pair.second);
                if (existingId != -1) {
                    if (m_pGraph.hasNode(existingId) &&
                        m_pGraph.getNode(existingId).type == PlanarizedGraph::NodeType::CROSSING) {
                        m_pGraph.destroyCrossing(existingId);
                    }
                }
                continue;
            }

            auto intersection = intersectOriginalEdgesForMove(pair.first, pair.second, -1, 0.0, 0.0);
            const int existingId = findExistingCrossingNodeForPair(pair.first, pair.second);

            if (intersection.has_value()) {
                if (existingId == -1) {
                    m_pGraph.createCrossing(pair.first, pair.second, intersection->first, intersection->second);
                } else {
                    if (m_pGraph.hasNode(existingId) &&
                        m_pGraph.getNode(existingId).type == PlanarizedGraph::NodeType::CROSSING) {
                        m_pGraph.updateNodePosition(existingId, intersection->first, intersection->second);
                    }
                }
            } else {
                if (existingId != -1) {
                    if (m_pGraph.hasNode(existingId) &&
                        m_pGraph.getNode(existingId).type == PlanarizedGraph::NodeType::CROSSING) {
                        m_pGraph.destroyCrossing(existingId);
                    }
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

    applyCrossingRemovals(plan);
    m_pGraph.updateNodePosition(variableNodeId, movedX, movedY);

    std::vector<int> crossingsToRefresh;
    crossingsToRefresh.reserve(32);
    m_pGraph.forEachNode([&](int nodeId, const PlanarizedGraph::PlanarNode& node) {
        if (node.type != PlanarizedGraph::NodeType::CROSSING) return;
        if (movedOriginalEdges.count(node.original_edge_1) == 0 &&
            movedOriginalEdges.count(node.original_edge_2) == 0) {
            return;
        }
        crossingsToRefresh.push_back(nodeId);
    });

    for (int crossingId : crossingsToRefresh) {
        if (!m_pGraph.hasNode(crossingId)) continue;
        const auto& crossing = m_pGraph.getNode(crossingId);
        auto newPos = intersectOriginalEdgesForMove(
            crossing.original_edge_1,
            crossing.original_edge_2,
            variableNodeId,
            movedX,
            movedY
        );
        if (!newPos.has_value()) continue;
        m_pGraph.updateNodePosition(crossingId, newPos->first, newPos->second);
    }

    applyCrossingInsertions(variableNodeId, movedX, movedY, plan);
    reconcileCrossingsForMovedEdges(movedOriginalEdges);

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
