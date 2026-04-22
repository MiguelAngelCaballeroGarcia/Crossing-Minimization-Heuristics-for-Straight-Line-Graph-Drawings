#pragma once

enum class LocalSegmentType {
    ORIGINAL_SUBSEGMENT,
    ROI_BORDER,
    RAY
};

/**
 * @brief A lightweight segment used exclusively for local region building.
 * It is guaranteed to be fully contained within the ROI.
 */
struct LocalSegment {
    double x1, y1;
    double x2, y2;

    // Links back to the originating planar/original edge.
    // Use -1 for ROI border walls and future emitted rays.
    int originalEdgeId;
    LocalSegmentType type;
    
    // For RAY type: the node ID that the ray is shot from.
    // For other types: -1 (not applicable).
    int raySourceNodeId = -1;
};