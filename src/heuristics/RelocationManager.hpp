#pragma once

#include <vector>
#include <optional>
#include <utility>
#include <random>
#include <unordered_set>
#include "../logic/PlanarizedGraph.hpp"
#include "LocalSegment.hpp"
#include "RegionBuilder.hpp"

/**
 * @brief Represents the physical and logical boundaries of our "frozen" grid cells.
 */
struct RegionOfInterest {
    // Grid indices for fast SpatialGrid lookups
    int minCol, maxCol;
    int minRow, maxRow;
    
    // Physical coordinates forming the absolute "Border Wall"
    double minX, maxX; 
    double minY, maxY;
};

/**
 * @brief Results from local region analysis including face weighting
 */
struct LocalRegionAnalysis {
    std::vector<LocalSegment> localGeometry;
    DualGraph dualGraph;
    std::vector<DualGraphEdge> dualTreeEdges;
    std::vector<int> faceParent;
    std::vector<double> faceWeights; // Weight for each face
    std::vector<int> faceGlobalCrossingDelta; // Exact global crossing delta for moving the node to each face
    int sourceFaceId = -1; // The face containing the selected node
    bool hasAmbiguousActiveSideProbe = false;
};

struct CrossingUpdatePlan {
    // Normalized original-edge pairs (min, max)
    std::vector<std::pair<int, int>> disappearingPairs;
    std::vector<std::pair<int, int>> appearingPairs;

    // Existing crossing nodes corresponding to disappearingPairs
    std::vector<int> crossingNodeIdsToDestroy;
};

struct RelocationStepResult {
    bool valid = false;
    bool moved = false;
    bool usedExplorationMove = false;
    int variableNodeId = -1;
    int sourceFaceId = -1;
    int targetFaceId = -1;
    double targetWeight = 0.0;
    int targetGlobalDelta = 0;
};

class RelocationManager {
public:
    // Takes references to the core logic structures
    RelocationManager(PlanarizedGraph& pGraph);

    /**
     * @brief Phase 1.1: Determines which ORIGINAL node should be relocated next.
     * @return The ID of the selected node.
     */
    int selectVariableNode();

    /**
     * @brief Phase 1.2: Finds the bounding box covering the node and its neighbors.
     * @param nodeId The node chosen for relocation.
     * @return The calculated ROI containing both grid indices and physical bounds.
     */
    RegionOfInterest calculateROI(int nodeId);

    /**
     * @brief Phase 1.3: Extracts all edges in the ROI, clips them to the borders,
     * and appends the ROI borders themselves as temporary segments.
     * @param roi The target region.
     * @param variableNodeId The node being moved (its incident edges might be ignored or handled specifically).
     * @return A flat, self-contained list of clipped segments ready for region subdivision.
     */
    std::vector<LocalSegment> extractAndClipGeometry(const RegionOfInterest& roi, int variableNodeId);

    /**
     * @brief Phase 1.4: Appends candidate rays for the variable node inside the ROI.
     * @param localGeometry Existing local segments (clipped edges + ROI walls).
     * @param roi The active region of interest.
     * @param variableNodeId The node being relocated.
     */
    void appendRaysToLocalGeometry(std::vector<LocalSegment>& localGeometry,
                                            const RegionOfInterest& roi,
                                            int variableNodeId);

    /**
     * @brief Analyzes local regions and computes face weights based on the weighting algorithm.
     * @param roi The region of interest.
     * @param variableNodeId The selected node for relocation.
     * @return Analysis with faces, dualgraph, and computed weights.
     */
    LocalRegionAnalysis analyzeLocalRegions(const RegionOfInterest& roi, int variableNodeId);

    /**
     * @brief Builds the set of crossing changes required to move the selected node
     * from its source face to targetFaceId by traversing the precomputed dual tree.
     */
    CrossingUpdatePlan buildCrossingUpdatePlan(int variableNodeId,
                                               const LocalRegionAnalysis& analysis,
                                               int targetFaceId) const;

    /**
     * @brief Removes crossings listed in the update plan.
     * Call this in the relocation transaction before inserting new crossings.
     */
    void applyCrossingRemovals(const CrossingUpdatePlan& plan);

    /**
     * @brief Inserts new crossings listed in the update plan.
     * Assumes the selected node position and incident original-edge geometry are already updated.
     */
    void applyCrossingInsertions(int variableNodeId,
                                 double movedX,
                                 double movedY,
                                 const CrossingUpdatePlan& plan);

    /**
     * @brief Executes one relocation update step.
     * Policy:
    * - Move to one of the minimum-weight faces.
     * - Ties are broken uniformly at random.
     */
    RelocationStepResult performRelocationStep();

private:
    PlanarizedGraph& m_pGraph;
    std::vector<int> m_originalNodeIds;
    std::unordered_set<int> m_originalNodeIdSet;  // For O(1) membership tests
    std::mt19937 m_rng;  // Seeded once in constructor for consistent performance

    /**
     * @brief Helper function to clip a segment against the ROI physical bounds.
     * Uses an algorithm like Cohen-Sutherland to trim the edge.
     * @return The clipped segment, or std::nullopt if the edge is completely outside.
     */
    std::optional<LocalSegment> clipToBoundingBox(double x1, double y1, double x2, double y2, 
                                                  const RegionOfInterest& roi, int edgeId);

    // Helpers for weighting algorithm
    int findSourceFace(const DualGraph& dualGraph, double x, double y) const;
    std::vector<double> computeFaceWeights(const DualGraph& dualGraph,
                                           int sourceFaceId,
                                           int variableNodeId,
                                           std::vector<DualGraphEdge>& treeEdges,
                                           std::vector<int>& faceParent,
                                           bool& ambiguousProbeDetected) const;

    std::pair<int, int> normalizeEdgePair(int edgeA, int edgeB) const;
    int computeActiveSideForBoundary(int currentFaceId,
                                     const DualGraphEdge& boundary,
                                     const DualGraph& dualGraph) const;
    std::vector<int> collectPathFromSource(int sourceFaceId,
                                           int targetFaceId,
                                           const std::vector<int>& faceParent) const;
    std::vector<DualGraphEdge> collectBoundaryPathEdges(const LocalRegionAnalysis& analysis,
                                                        int targetFaceId) const;

    std::vector<std::pair<int, int>> collectSelectedNeighborEdges(int variableNodeId) const;
    std::vector<std::pair<int, int>> collectIncidentEdgesForNode(int nodeId) const;

    void collectTransitionCrossingPairs(int variableNodeId,
                                        int currentFaceId,
                                        const DualGraphEdge& boundary,
                                        std::vector<std::pair<int, int>>& disappear,
                                        std::vector<std::pair<int, int>>& appear,
                                        const DualGraph& dualGraph) const;

    std::optional<std::pair<int, int>> findOriginalEdgeEndpoints(int originalEdgeId) const;
    std::optional<std::pair<double, double>> intersectOriginalEdgesForMove(int edgeA,
                                                                           int edgeB,
                                                                           int variableNodeId,
                                                                           double movedX,
                                                                           double movedY) const;
    int findExistingCrossingNodeForPair(int edgeA, int edgeB) const;
    int countGlobalCrossingsForVariableAtPosition(int variableNodeId,
                                                  double movedX,
                                                  double movedY,
                                                  const std::unordered_set<int>& variableIncidentEdges) const;
    int evaluateGlobalCrossingDeltaForMove(int variableNodeId,
                                           double movedX,
                                           double movedY,
                                           const std::unordered_set<int>& variableIncidentEdges) const;

    std::optional<std::pair<int, double>> chooseTargetFace(const LocalRegionAnalysis& analysis,
                                                           std::mt19937& rng) const;
    std::optional<std::pair<double, double>> chooseInteriorPointInFace(const Face& face) const;
    void reconcileCrossingsForMovedEdges(const std::unordered_set<int>& movedOriginalEdges);
    std::vector<int> computeFaceGlobalCrossingDeltas(const DualGraph& dualGraph,
                                                     int variableNodeId) const;
    
    // Determine which side of an infinite directed line p1->p2 a point falls on.
    int whichSideOfLine(double px, double py, double x1, double y1, double x2, double y2) const;
    
};