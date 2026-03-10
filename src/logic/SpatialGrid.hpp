#pragma once

#include <vector>
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
    SpatialGrid(double min_x, double max_x, double min_y, double max_y, int totalNodes);

    void insertNode(int nodeIdx, double x, double y);
    
    // Remember: an edge might be inserted into multiple cells!
    void insertEdge(int edgeIdx, double x1, double y1, double x2, double y2);

    // Get a cell for querying
    const GridCell& getCell(int index) const { return cells[index]; }
};