#pragma once
#include <vector>
#include "../graph/graph.hpp"

class IntersectionDetector {
public:
    IntersectionDetector() = default;
    int countIntersections(const Graph& G);
};
