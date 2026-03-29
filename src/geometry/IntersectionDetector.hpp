#pragma once
#include <vector>
#include "../logic/Graph.hpp"
#include "../logic/SpatialGrid.hpp"
#include "../logic/PlanarizedGraph.hpp"

class IntersectionDetector {
public:
	static std::vector<PlanarizedGraph::IntersectionData> findIntersections(const Graph& graph, const SpatialGrid& grid);
};
