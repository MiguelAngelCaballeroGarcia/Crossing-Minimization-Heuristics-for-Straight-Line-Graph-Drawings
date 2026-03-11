#pragma once

#include <vector>
#include <algorithm>
#include "Graph.hpp"

struct GridCell {
    std::vector<int> nodeIndices; // Indices into Graph::nodes
    std::vector<int> edgeIndices; // Indices into Graph::edges
};

class SpatialGrid {
private:
    double minX, maxX, minY, maxY;
    double cellWidth, cellHeight;
    int numCellsX, numCellsY;
    
    // The actual storage: a 1D vector representing the 2D grid
    std::vector<GridCell> cells;

    // Helper to convert (x, y) coordinates to a 1D grid index
    int getCellIndex(double x, double y) const;

public:
// Bulk-loading constructor
    SpatialGrid(const Graph& graph);

    // Getters for the PlanarizedGraph to build the Floating Frame
    double getMinX() const { return minX; };
    double getMaxX() const { return maxX; };
    double getMinY() const { return minY; };
    double getMaxY() const { return maxY; };

    // Public Incremental Update Methods
    void insertNode(int nodeIdx, double x, double y);
    void insertEdge(int edgeIdx, double x1, double y1, double x2, double y2);
    
    // NEW: Required for incremental updates when nodes move
    void removeEdge(int edgeIdx, double x1, double y1, double x2, double y2);

    const GridCell& getCell(int index) const { return cells[index]; };
};