#include "RegionBuilder.hpp"
#include <cmath>
#include <algorithm>
#include <functional>

namespace {
    struct SegmentBounds {
        double minX;
        double maxX;
        double minY;
        double maxY;
    };

    SegmentBounds getBounds(const LocalSegment& segment) {
        return {
            std::min(segment.x1, segment.x2),
            std::max(segment.x1, segment.x2),
            std::min(segment.y1, segment.y2),
            std::max(segment.y1, segment.y2)
        };
    }

    bool boundsOverlap(const SegmentBounds& a, const SegmentBounds& b) {
        return !(a.maxX < b.minX || b.maxX < a.minX || a.maxY < b.minY || b.maxY < a.minY);
    }
}

std::size_t RegionBuilder::PointHash::operator()(const Point& p) const noexcept {
    const long long quantizedX = static_cast<long long>(std::llround(p.x / EPSILON));
    const long long quantizedY = static_cast<long long>(std::llround(p.y / EPSILON));
    std::size_t seed = 0;
    seed ^= std::hash<long long>{}(quantizedX) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<long long>{}(quantizedY) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

// Helper to get intersection of two segments
static bool getIntersection(double x1, double y1, double x2, double y2,
                            double x3, double y3, double x4, double y4,
                            double& ix, double& iy) {
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (std::abs(denom) < RegionBuilder::EPSILON) return false;

    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
    double u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;

    if (t >= -RegionBuilder::EPSILON && t <= 1.0 + RegionBuilder::EPSILON &&
        u >= -RegionBuilder::EPSILON && u <= 1.0 + RegionBuilder::EPSILON) {
        ix = x1 + t * (x2 - x1);
        iy = y1 + t * (y2 - y1);
        return true;
    }
    return false;
}

std::vector<LocalSegment> RegionBuilder::splitSegmentsAtIntersections(const std::vector<LocalSegment>& segments) {
    std::vector<LocalSegment> result;

    // PRECOMPUTE BOUNDS HERE
    std::vector<SegmentBounds> bounds(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
        bounds[i] = getBounds(segments[i]);
    }

    // O(N^2) intersection finding. Safe because N is bounded by the small local ROI.
    for (size_t i = 0; i < segments.size(); ++i) {
        const auto boundsI = getBounds(segments[i]);
        std::vector<Point> splits;
        splits.push_back({segments[i].x1, segments[i].y1});
        splits.push_back({segments[i].x2, segments[i].y2});

        for (size_t j = 0; j < segments.size(); ++j) {
            if (i == j) continue;
            if (!boundsOverlap(bounds[i], bounds[j])) continue;
            double ix, iy;
            if (getIntersection(segments[i].x1, segments[i].y1, segments[i].x2, segments[i].y2,
                                segments[j].x1, segments[j].y1, segments[j].x2, segments[j].y2, ix, iy)) {
                splits.push_back({ix, iy});
            } else {
                // Collinearity fallback: Check if endpoints of J lie on segment I
                auto isPointOnSegment = [](double px, double py, const LocalSegment& seg) {
                    double cross = (px - seg.x1) * (seg.y2 - seg.y1) - (py - seg.y1) * (seg.x2 - seg.x1);
                    if (std::abs(cross) > RegionBuilder::EPSILON) return false; // Not collinear
                    
                    // Check if within bounds
                    double minX = std::min(seg.x1, seg.x2), maxX = std::max(seg.x1, seg.x2);
                    double minY = std::min(seg.y1, seg.y2), maxY = std::max(seg.y1, seg.y2);
                    return px >= minX - RegionBuilder::EPSILON && px <= maxX + RegionBuilder::EPSILON &&
                        py >= minY - RegionBuilder::EPSILON && py <= maxY + RegionBuilder::EPSILON;
                };

                if (isPointOnSegment(segments[j].x1, segments[j].y1, segments[i])) splits.push_back({segments[j].x1, segments[j].y1});
                if (isPointOnSegment(segments[j].x2, segments[j].y2, segments[i])) splits.push_back({segments[j].x2, segments[j].y2});
            }
        }

        // Sort points along the segment to create ordered sub-segments.
        // Using the segment direction keeps the ordering stable for collinear points.
        const Point p1 = splits[0];
        const Point p2 = splits[1];
        const double dx = p2.x - p1.x;
        const double dy = p2.y - p1.y;
        const double len2 = dx * dx + dy * dy;

        if (len2 < EPSILON * EPSILON) {
            continue;
        }

        std::sort(splits.begin(), splits.end(), [&p1, dx, dy, len2](const Point& a, const Point& b) {
            const double ta = ((a.x - p1.x) * dx + (a.y - p1.y) * dy) / len2;
            const double tb = ((b.x - p1.x) * dx + (b.y - p1.y) * dy) / len2;
            return ta < tb;
        });

        splits.erase(std::unique(splits.begin(), splits.end(), [](const Point& a, const Point& b) {
            return a == b;
        }), splits.end());

        for (size_t k = 0; k < splits.size() - 1; ++k) {
            if (!(splits[k] == splits[k+1])) {
                LocalSegment subseg;
                subseg.x1 = splits[k].x;
                subseg.y1 = splits[k].y;
                subseg.x2 = splits[k + 1].x;
                subseg.y2 = splits[k + 1].y;
                subseg.originalEdgeId = segments[i].originalEdgeId;
                subseg.type = segments[i].type;
                subseg.raySourceNodeId = segments[i].raySourceNodeId;
                result.push_back(subseg);
            }
        }
    }
    return result;
}

DualGraph RegionBuilder::buildRegionsAndDualGraph(const std::vector<LocalSegment>& inputSegments) {
    auto segments = splitSegmentsAtIntersections(inputSegments);

    std::vector<Point> vertices;
    std::unordered_map<Point, int, PointHash> vertexLookup;
    std::vector<HalfEdge> halfEdges;

    // 1. Extract unique vertices and create twins
    auto getVertexId = [&vertices, &vertexLookup](double x, double y) {
        Point p{x, y};
        auto found = vertexLookup.find(p);
        if (found != vertexLookup.end()) {
            return found->second;
        }
        vertices.push_back(p);
        int newId = static_cast<int>(vertices.size() - 1);
        vertexLookup.emplace(p, newId);
        return newId;
    };

    for (const auto& seg : segments) {
        int u = getVertexId(seg.x1, seg.y1);
        int v = getVertexId(seg.x2, seg.y2);

        int he1_id = halfEdges.size();
        int he2_id = halfEdges.size() + 1;

        halfEdges.push_back({he1_id, u, he2_id, -1, -1, seg.originalEdgeId, seg.type, seg.raySourceNodeId});
        halfEdges.push_back({he2_id, v, he1_id, -1, -1, seg.originalEdgeId, seg.type, seg.raySourceNodeId});
    }

    // 2. Link Next pointers via angular sorting
    // Group outgoing half-edges by their origin vertex
    std::vector<std::vector<int>> vertexOutgoing(vertices.size());
    for (const auto& he : halfEdges) {
        vertexOutgoing[he.originVertexId].push_back(he.id);
    }

    for (size_t i = 0; i < vertices.size(); ++i) {
        auto& outEdges = vertexOutgoing[i];
        const auto& origin = vertices[i];

        // Sort radially based on angle of destination vertex
        std::sort(outEdges.begin(), outEdges.end(), [&](int aId, int bId) {
            int destA = halfEdges[halfEdges[aId].twinId].originVertexId;
            int destB = halfEdges[halfEdges[bId].twinId].originVertexId;
            double angleA = std::atan2(vertices[destA].y - origin.y, vertices[destA].x - origin.x);
            double angleB = std::atan2(vertices[destB].y - origin.y, vertices[destB].x - origin.x);
            return angleA < angleB;
        });

        // Link twins to next outgoing
        for (size_t e = 0; e < outEdges.size(); ++e) {
            int currentOut = outEdges[e];
            int nextOut = outEdges[(e + 1) % outEdges.size()];
            // The twin of current comes INTO the vertex. Its 'next' is the next outgoing edge counter-clockwise.
            halfEdges[halfEdges[currentOut].twinId].nextId = nextOut;
        }
    }

    // 3. Extract Faces
    DualGraph dualGraph;
    int currentFaceId = 0;
    std::vector<double> signedAreas;

    for (auto& he : halfEdges) {
        if (he.faceId != -1) continue; // Already visited

        Face face;
        face.id = currentFaceId++;
        int startId = he.id;
        int currId = startId;
        double area = 0.0;

        do {
            halfEdges[currId].faceId = face.id;
            face.halfEdgeIds.push_back(currId);

            // Collect face vertex
            Point p1 = vertices[halfEdges[currId].originVertexId];
            face.vertices.push_back({p1.x, p1.y});

            // Calculate polygon area (Shoelace formula) to find outer face
            Point p2 = vertices[halfEdges[halfEdges[currId].twinId].originVertexId];
            area += (p1.x * p2.y - p2.x * p1.y);

            currId = halfEdges[currId].nextId;
        } while (currId != startId);

        // Temporary value. Outer face is assigned globally below.
        face.isOuter = false;
        signedAreas.push_back(area);

        dualGraph.faces.push_back(face);
    }

    // Assign a single outer face: the one with largest absolute area.
    int outerFaceId = -1;
    double maxAbsArea = -1.0;
    for (size_t i = 0; i < signedAreas.size(); ++i) {
        const double absArea = std::abs(signedAreas[i]);
        if (absArea > maxAbsArea) {
            maxAbsArea = absArea;
            outerFaceId = static_cast<int>(i);
        }
    }

    if (outerFaceId >= 0 && outerFaceId < static_cast<int>(dualGraph.faces.size())) {
        dualGraph.faces[outerFaceId].isOuter = true;
    }

    // 4. Construct Dual Graph Adjacency
    for (const auto& face : dualGraph.faces) {
        if (face.isOuter) continue;

        for (int heId : face.halfEdgeIds) {
            int twinFaceId = halfEdges[halfEdges[heId].twinId].faceId;

            // Only create an edge if it connects two different inner faces
            if (face.id < twinFaceId && !dualGraph.faces[twinFaceId].isOuter) {
                const auto& p1 = vertices[halfEdges[heId].originVertexId];
                const auto& p2 = vertices[halfEdges[halfEdges[heId].twinId].originVertexId];
                dualGraph.adjacency.push_back({
                    face.id,
                    twinFaceId,
                    halfEdges[heId].originalEdgeId,
                    halfEdges[heId].boundaryType,
                    halfEdges[heId].raySourceNodeId,
                    p1.x,
                    p1.y,
                    p2.x,
                    p2.y
                });
            }
        }
    }

    return dualGraph;
}