#pragma once

#include <vector>
#include <algorithm>

struct GridCell {
    std::vector<int> nodeIndices; 
    std::vector<int> edgeIndices; 
};

class SpatialGrid {
private:
    double minX, maxX, minY, maxY;
    double cellWidth, cellHeight;
    int numCellsX, numCellsY;
    
    std::vector<GridCell> cells;

    // Helper to convert (x, y) coordinates to a 1D grid index
    int getCellIndex(double x, double y) const;

public:
    // Default constructor (required so PlanarizedGraph can declare it as a member and initialize it later)
    SpatialGrid();

    // Generic Geometric Constructor
    SpatialGrid(double minX, double maxX, double minY, double maxY, int expectedTotalNodes);

    // Getters for the Floating Frame
    double getMinX() const { return minX; }
    double getMaxX() const { return maxX; }
    double getMinY() const { return minY; }
    double getMaxY() const { return maxY; }

    size_t getNumCells() const { return cells.size(); }

    double getCellWidth() const { return cellWidth; }
    double getCellHeight() const { return cellHeight; }
    int getNumCellsX() const { return numCellsX; }
    int getNumCellsY() const { return numCellsY; }

    // Public Incremental Update Methods
    void insertNode(int nodeIdx, double x, double y);
    void removeNode(int nodeIdx, double x, double y);

    void insertEdge(int edgeIdx, double x1, double y1, double x2, double y2);
    void removeEdge(int edgeIdx, double x1, double y1, double x2, double y2);

    const GridCell& getCell(int index) const { return cells[index]; }
};