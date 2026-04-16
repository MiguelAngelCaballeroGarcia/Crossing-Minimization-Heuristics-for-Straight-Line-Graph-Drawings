#pragma once

#include <vector>
#include <unordered_set>

#include "../RelocationManager.hpp"

namespace relocation_detail {

int computeOutCode(double x, double y, const RegionOfInterest& roi);
int toColIndex(double x, const SpatialGrid& grid);
int toRowIndex(double y, const SpatialGrid& grid);
int toCellIndex(int col, int row, const SpatialGrid& grid);

std::vector<int> collectOriginalNeighbors(int nodeId, const PlanarizedGraph& pGraph);
std::unordered_set<int> collectIncidentOriginalEdgeIds(int nodeId, const PlanarizedGraph& pGraph);

bool pointInPolygon(double x, double y, const std::vector<FaceVertex>& polygon);

} // namespace relocation_detail
