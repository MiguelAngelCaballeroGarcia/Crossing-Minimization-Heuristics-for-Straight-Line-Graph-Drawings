#pragma once

#include <vector>
#include <cmath>
#include <cstddef>
#include <unordered_map>
#include "LocalSegment.hpp"

struct FaceVertex {
    double x;
    double y;
};

struct DualGraphEdge {
    int faceA;
    int faceB;
    int originalEdgeId; // The physical edge separating them (-1 for rays/borders)
    LocalSegmentType boundaryType; // ORIGINAL_SUBSEGMENT or RAY
    int raySourceNodeId; // For RAY type: which node the ray came from; otherwise -1
    int rayEmittingNodeId; // For RAY type: which node the ray is emitted from; otherwise -1
    double x1;
    double y1;
    double x2;
    double y2;
};

struct Face {
    int id;
    std::vector<int> halfEdgeIds;
    std::vector<FaceVertex> vertices; // Vertices of the face polygon
    bool isOuter;
};

struct DualGraph {
    std::vector<Face> faces;
    std::vector<DualGraphEdge> adjacency;
};

class RegionBuilder {
public:

    static constexpr double EPSILON = 1e-7;

    /**
     * @brief Takes raw segments, splits them at intersections, builds a DCEL,
     * extracts faces, and constructs the dual graph.
     */
    static DualGraph buildRegionsAndDualGraph(const std::vector<LocalSegment>& segments);

private:

    struct Point {
        double x, y;
        bool operator==(const Point& o) const { return std::abs(x - o.x) < EPSILON && std::abs(y - o.y) < EPSILON; }

    };

    struct PointHash {
        std::size_t operator()(const Point& p) const noexcept;
    };

    struct HalfEdge {
        int id;
        int originVertexId;
        int twinId = -1;
        int nextId = -1;
        int faceId = -1;
        int originalEdgeId;
        LocalSegmentType boundaryType; // Store segment type
        int raySourceNodeId = -1; // For RAY type: which node it came from
        int rayEmittingNodeId = -1; // For RAY type: which node determines the ray direction
    };

    static std::vector<LocalSegment> splitSegmentsAtIntersections(const std::vector<LocalSegment>& segments);
};