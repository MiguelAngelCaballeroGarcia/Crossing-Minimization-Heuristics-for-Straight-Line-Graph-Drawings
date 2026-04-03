#pragma once

#include <vector>
#include <optional>
#include "../logic/PlanarizedGraph.hpp"

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
 * @brief A lightweight segment used exclusively for local region building.
 * It is guaranteed to be fully contained within the ROI.
 */
struct LocalSegment {
    double x1, y1;
    double x2, y2;
    
    // Links back to the PlanarizedGraph. 
    // Use -1 for ROI border walls and future emitted rays.
    int originalEdgeId; 
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

private:
    PlanarizedGraph& m_pGraph;
    std::vector<int> m_originalNodeIds;

    /**
     * @brief Helper function to clip a segment against the ROI physical bounds.
     * Uses an algorithm like Cohen-Sutherland to trim the edge.
     * @return The clipped segment, or std::nullopt if the edge is completely outside.
     */
    std::optional<LocalSegment> clipToBoundingBox(double x1, double y1, double x2, double y2, 
                                                  const RegionOfInterest& roi, int edgeId);
};