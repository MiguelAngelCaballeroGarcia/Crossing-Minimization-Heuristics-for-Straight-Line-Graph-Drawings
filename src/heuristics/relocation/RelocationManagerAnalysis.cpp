#include "../RelocationManager.hpp"
#include "RelocationHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>

LocalRegionAnalysis RelocationManager::analyzeLocalRegions(const RegionOfInterest& roi, int variableNodeId) {
    LocalRegionAnalysis analysis;
    resetStepCrossingCaches();

    analysis.localGeometry = extractAndClipGeometry(roi, variableNodeId);
    analysis.dualGraph = RegionBuilder::buildRegionsAndDualGraph(analysis.localGeometry);

    if (!m_pGraph.hasNode(variableNodeId)) {
        analysis.sourceFaceId = -1;
        analysis.faceWeights.resize(analysis.dualGraph.faces.size(), 0.0);
        return analysis;
    }

    const auto& variableNode = m_pGraph.getNode(variableNodeId);
    analysis.sourceFaceId = findSourceFace(analysis.dualGraph, variableNode.x, variableNode.y);

    bool ambiguousProbeDetected = false;
    analysis.faceWeights = computeFaceWeights(analysis.dualGraph,
                                              analysis.sourceFaceId,
                                              variableNodeId,
                                              analysis.dualTreeEdges,
                                              analysis.faceParent,
                                              ambiguousProbeDetected);
    analysis.hasAmbiguousActiveSideProbe = ambiguousProbeDetected;
    analysis.faceGlobalCrossingDelta = computeFaceGlobalCrossingDeltas(analysis, variableNodeId);

    return analysis;
}

std::vector<int> RelocationManager::computeFaceGlobalCrossingDeltas(const LocalRegionAnalysis& analysis,
                                                                    int variableNodeId) const {
    std::vector<int> deltas(analysis.dualGraph.faces.size(), 0);

    if (analysis.sourceFaceId < 0 ||
        analysis.sourceFaceId >= static_cast<int>(analysis.dualGraph.faces.size())) {
        return deltas;
    }

    auto makeFacePairKey = [](int a, int b) {
        const std::uint64_t lo = static_cast<std::uint64_t>(std::min(a, b));
        const std::uint64_t hi = static_cast<std::uint64_t>(std::max(a, b));
        return (lo << 32) | hi;
    };

    std::unordered_map<std::uint64_t, DualGraphEdge> treeEdgeByFaces;
    treeEdgeByFaces.reserve(analysis.dualTreeEdges.size() * 2);
    for (const auto& edge : analysis.dualTreeEdges) {
        treeEdgeByFaces[makeFacePairKey(edge.faceA, edge.faceB)] = edge;
    }

    auto findBoundaryEdge = [&](int faceA, int faceB) -> std::optional<DualGraphEdge> {
        const auto key = makeFacePairKey(faceA, faceB);
        auto treeIt = treeEdgeByFaces.find(key);
        if (treeIt != treeEdgeByFaces.end()) {
            return treeIt->second;
        }

        for (const auto& edge : analysis.dualGraph.adjacency) {
            if ((edge.faceA == faceA && edge.faceB == faceB) ||
                (edge.faceA == faceB && edge.faceB == faceA)) {
                return edge;
            }
        }
        return std::nullopt;
    };

    std::vector<std::vector<int>> children(analysis.faceParent.size());
    for (int faceId = 0; faceId < static_cast<int>(analysis.faceParent.size()); ++faceId) {
        const int parent = analysis.faceParent[faceId];
        if (parent >= 0 && parent < static_cast<int>(children.size())) {
            children[parent].push_back(faceId);
        }
    }

    std::queue<int> q;
    std::vector<std::uint8_t> visited(analysis.dualGraph.faces.size(), 0);
    q.push(analysis.sourceFaceId);
    visited[analysis.sourceFaceId] = 1;

    while (!q.empty()) {
        const int currentFace = q.front();
        q.pop();

        if (currentFace < 0 || currentFace >= static_cast<int>(children.size())) continue;

        for (int childFace : children[currentFace]) {
            auto boundary = findBoundaryEdge(currentFace, childFace);
            if (!boundary.has_value()) {
                continue;
            }

            std::vector<std::pair<int, int>> disappear;
            std::vector<std::pair<int, int>> appear;
            collectTransitionCrossingPairs(variableNodeId,
                                           currentFace,
                                           *boundary,
                                           disappear,
                                           appear,
                                           analysis.dualGraph);

            std::set<std::pair<int, int>> disappearSet;
            std::set<std::pair<int, int>> appearSet;
            for (const auto& p : disappear) disappearSet.insert(normalizeEdgePair(p.first, p.second));
            for (const auto& p : appear) appearSet.insert(normalizeEdgePair(p.first, p.second));

            int boundaryDelta = 0;
            for (const auto& p : appearSet) {
                if (disappearSet.find(p) == disappearSet.end()) {
                    ++boundaryDelta;
                }
            }
            for (const auto& p : disappearSet) {
                if (appearSet.find(p) == appearSet.end()) {
                    --boundaryDelta;
                }
            }

            deltas[childFace] = deltas[currentFace] + boundaryDelta;
            if (!visited[childFace]) {
                visited[childFace] = 1;
                q.push(childFace);
            }
        }
    }

    return deltas;
}

int RelocationManager::findSourceFace(const DualGraph& dualGraph, double x, double y) const {
    for (const auto& face : dualGraph.faces) {
        if (face.isOuter) continue;
        if (face.vertices.empty()) continue;

        if (relocation_detail::pointInPolygon(x, y, face.vertices)) {
            return face.id;
        }
    }

    for (const auto& face : dualGraph.faces) {
        if (!face.isOuter) {
            return face.id;
        }
    }

    return -1;
}

int RelocationManager::whichSideOfLine(double px, double py, double x1, double y1, double x2, double y2) const {
    const double dx = x2 - x1;
    const double dy = y2 - y1;
    const double toPx = px - x1;
    const double toPy = py - y1;
    const double cross = dx * toPy - dy * toPx;
    if (std::abs(cross) < 1e-11) return 0;
    return (cross > 0.0) ? 1 : -1;
}

std::vector<double> RelocationManager::computeFaceWeights(const DualGraph& dualGraph,
                                                          int sourceFaceId,
                                                          int variableNodeId,
                                                          std::vector<DualGraphEdge>& treeEdges,
                                                          std::vector<int>& faceParent,
                                                          bool& ambiguousProbeDetected) const {
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> weights(dualGraph.faces.size(), INF);
    faceParent.assign(dualGraph.faces.size(), -1);
    treeEdges.clear();

    if (sourceFaceId < 0 || sourceFaceId >= static_cast<int>(dualGraph.faces.size())) {
        return weights;
    }

    auto variableNeighbors = relocation_detail::collectOriginalNeighbors(variableNodeId, m_pGraph);

    std::vector<std::vector<int>> faceAdjEdges(dualGraph.faces.size());
    for (size_t edgeIdx = 0; edgeIdx < dualGraph.adjacency.size(); ++edgeIdx) {
        const auto& e = dualGraph.adjacency[edgeIdx];
        if (e.faceA >= 0 && e.faceA < static_cast<int>(faceAdjEdges.size())) {
            faceAdjEdges[e.faceA].push_back(static_cast<int>(edgeIdx));
        }
        if (e.faceB >= 0 && e.faceB < static_cast<int>(faceAdjEdges.size())) {
            faceAdjEdges[e.faceB].push_back(static_cast<int>(edgeIdx));
        }
    }

    std::queue<int> bfsQueue;
    std::vector<bool> visited(dualGraph.faces.size(), false);

    bfsQueue.push(sourceFaceId);
    weights[sourceFaceId] = 0.0;
    visited[sourceFaceId] = true;

    while (!bfsQueue.empty()) {
        int currentFace = bfsQueue.front();
        bfsQueue.pop();

        double currentWeight = weights[currentFace];

        for (int edgeIdx : faceAdjEdges[currentFace]) {
            const auto& dualEdge = dualGraph.adjacency[edgeIdx];
            int nextFace = -1;

            if (dualEdge.faceA == currentFace) {
                nextFace = dualEdge.faceB;
            } else if (dualEdge.faceB == currentFace) {
                nextFace = dualEdge.faceA;
            } else {
                continue;
            }

            if (visited[nextFace]) continue;

            double delta = 0.0;
            const auto& currentFacePolygon = dualGraph.faces[currentFace].vertices;

            const double dx = dualEdge.x2 - dualEdge.x1;
            const double dy = dualEdge.y2 - dualEdge.y1;
            const double len = std::sqrt(dx * dx + dy * dy);
            const double midX = (dualEdge.x1 + dualEdge.x2) * 0.5;
            const double midY = (dualEdge.y1 + dualEdge.y2) * 0.5;
            const double perpX = (len < 1e-11) ? 0.0 : -dy / len;
            const double perpY = (len < 1e-11) ? 0.0 : dx / len;
            const double epsilon = 1e-7;

            const std::pair<double, double> candidateA{midX + epsilon * perpX, midY + epsilon * perpY};
            const std::pair<double, double> candidateB{midX - epsilon * perpX, midY - epsilon * perpY};

            const bool aInsideCurrentFace = relocation_detail::pointInPolygon(
                candidateA.first,
                candidateA.second,
                currentFacePolygon
            );
            const bool bInsideCurrentFace = relocation_detail::pointInPolygon(
                candidateB.first,
                candidateB.second,
                currentFacePolygon
            );

            std::pair<double, double> activePoint = candidateA;
            if (!aInsideCurrentFace && bInsideCurrentFace) {
                activePoint = candidateB;
            } else if (!bInsideCurrentFace && aInsideCurrentFace) {
                activePoint = candidateA;
            } else {
                ambiguousProbeDetected = true;
            }

            int activeSide = whichSideOfLine(
                activePoint.first,
                activePoint.second,
                dualEdge.x1,
                dualEdge.y1,
                dualEdge.x2,
                dualEdge.y2
            );

            if (activeSide == 0) activeSide = 1;

            if (dualEdge.boundaryType == LocalSegmentType::ORIGINAL_SUBSEGMENT) {
                int sameSideCount = 0;
                int otherSideCount = 0;

                for (int neighborId : variableNeighbors) {
                    if (!m_pGraph.hasNode(neighborId)) continue;
                    const auto& neighborNode = m_pGraph.getNode(neighborId);

                    int side = whichSideOfLine(
                        neighborNode.x,
                        neighborNode.y,
                        dualEdge.x1,
                        dualEdge.y1,
                        dualEdge.x2,
                        dualEdge.y2
                    );
                    if (side == 0) continue;
                    if (side == activeSide) sameSideCount++;
                    else otherSideCount++;
                }

                delta = static_cast<double>(sameSideCount - otherSideCount);
            } else if (dualEdge.boundaryType == LocalSegmentType::RAY) {
                int raySourceId = dualEdge.raySourceNodeId;

                if (m_pGraph.hasNode(raySourceId)) {
                    auto raySourceNeighbors = relocation_detail::collectOriginalNeighbors(raySourceId, m_pGraph);

                    int sameSideCount = 0;
                    int otherSideCount = 0;

                    for (int neighborId : raySourceNeighbors) {
                        if (neighborId == variableNodeId) continue;
                        if (!m_pGraph.hasNode(neighborId)) continue;
                        const auto& neighborNode = m_pGraph.getNode(neighborId);

                        int side = whichSideOfLine(
                            neighborNode.x,
                            neighborNode.y,
                            dualEdge.x1,
                            dualEdge.y1,
                            dualEdge.x2,
                            dualEdge.y2
                        );

                        if (side == 0) continue;
                        if (side == activeSide) sameSideCount++;
                        else otherSideCount++;
                    }

                    delta = static_cast<double>(otherSideCount - sameSideCount);
                }
            }

            double newWeight = currentWeight + delta;
            weights[nextFace] = newWeight;
            visited[nextFace] = true;
            faceParent[nextFace] = currentFace;
            treeEdges.push_back(dualEdge);
            bfsQueue.push(nextFace);
        }
    }

    return weights;
}

std::optional<std::pair<int, double>> RelocationManager::chooseTargetFace(const LocalRegionAnalysis& analysis,
                                                                           std::mt19937& rng) const {
    std::vector<int> candidateFaces;
    candidateFaces.reserve(analysis.dualGraph.faces.size());

    for (const auto& face : analysis.dualGraph.faces) {
        if (face.isOuter) continue;
        if (face.id == analysis.sourceFaceId) continue;
        if (face.id < 0 || face.id >= static_cast<int>(analysis.faceWeights.size())) continue;
        if (!std::isfinite(analysis.faceWeights[face.id])) continue;
        candidateFaces.push_back(face.id);
    }

    if (candidateFaces.empty()) return std::nullopt;

    double bestDelta = std::numeric_limits<double>::infinity();
    for (int faceId : candidateFaces) {
        bestDelta = std::min(bestDelta, analysis.faceWeights[faceId]);
    }

    if (bestDelta > 0) {
        return std::nullopt;
    }

    std::vector<int> ties;
    for (int faceId : candidateFaces) {
        if (analysis.faceWeights[faceId] == bestDelta) {
            ties.push_back(faceId);
        }
    }

    std::uniform_int_distribution<size_t> pickTie(0, ties.size() - 1);
    const int chosenFace = ties[pickTie(rng)];
    return std::make_pair(chosenFace, bestDelta);
}

std::optional<std::pair<double, double>> RelocationManager::chooseInteriorPointInFace(const Face& face) const {
    if (face.vertices.empty()) return std::nullopt;

    double cx = 0.0;
    double cy = 0.0;
    for (const auto& v : face.vertices) {
        cx += v.x;
        cy += v.y;
    }
    cx /= static_cast<double>(face.vertices.size());
    cy /= static_cast<double>(face.vertices.size());

    return std::make_pair(cx, cy);
}
