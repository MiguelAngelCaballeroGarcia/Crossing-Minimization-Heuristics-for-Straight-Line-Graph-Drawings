#include "SpatialGrid.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

SpatialGrid::SpatialGrid() 
    : minX(0), maxX(0), minY(0), maxY(0), 
      cellWidth(1), cellHeight(1), numCellsX(0), numCellsY(0) {}

SpatialGrid::SpatialGrid(double minX, double maxX, double minY, double maxY, int expectedTotalNodes) 
    : minX(minX), maxX(maxX), minY(minY), maxY(maxY) 
{
    double width = maxX - minX;
    double height = maxY - minY;
    
    // Prevent divide-by-zero on perfectly flat graphs
    if (width == 0.0) width = 1.0;
    if (height == 0.0) height = 1.0;

    double targetTotalCells = static_cast<double>(expectedTotalNodes) / 5.0;
    if (targetTotalCells < 1.0) targetTotalCells = 1.0;

    double aspectRatio = width / height;
    numCellsX = static_cast<int>(std::sqrt(targetTotalCells * aspectRatio));
    numCellsX = std::max(1, numCellsX);
    
    numCellsY = static_cast<int>(std::ceil(targetTotalCells / numCellsX));
    numCellsY = std::max(1, numCellsY);

    cellWidth = width / numCellsX;
    cellHeight = height / numCellsY;

    cells.resize(numCellsX * numCellsY);
}

int SpatialGrid::getCellIndex(double x, double y) const {
    int ix = static_cast<int>(std::floor((x - minX) / cellWidth));
    int iy = static_cast<int>(std::floor((y - minY) / cellHeight));

    if (ix < 0) ix = 0;
    if (ix >= numCellsX) ix = numCellsX - 1;
    if (iy < 0) iy = 0;
    if (iy >= numCellsY) iy = numCellsY - 1;

    return iy * numCellsX + ix;
}

void SpatialGrid::insertNode(int nodeIdx, double x, double y) {
    if (cells.empty()) return;
    int cellIndex = getCellIndex(x, y);
    cells[cellIndex].nodeIndices.push_back(nodeIdx);
}

void SpatialGrid::removeNode(int nodeIdx, double x, double y) {
    if (cells.empty()) return;
    int cellIndex = getCellIndex(x, y);
    auto& nodeList = cells[cellIndex].nodeIndices;
    auto it = std::find(nodeList.begin(), nodeList.end(), nodeIdx);
    if (it != nodeList.end()) {
        nodeList.erase(it);
    }
}

void SpatialGrid::insertEdge(int edgeIdx, double x1, double y1, double x2, double y2) {
    if (cells.empty()) return;
    // --- SAME DDA LOGIC YOU WROTE ---
    int ix = static_cast<int>(std::floor((x1 - minX) / cellWidth));
    int iy = static_cast<int>(std::floor((y1 - minY) / cellHeight));
    int iendX = static_cast<int>(std::floor((x2 - minX) / cellWidth));
    int iendY = static_cast<int>(std::floor((y2 - minY) / cellHeight));

    ix = std::clamp(ix, 0, numCellsX - 1);
    iy = std::clamp(iy, 0, numCellsY - 1);
    iendX = std::clamp(iendX, 0, numCellsX - 1);
    iendY = std::clamp(iendY, 0, numCellsY - 1);

    int stepX = (x2 > x1) ? 1 : (x2 < x1) ? -1 : 0;
    int stepY = (y2 > y1) ? 1 : (y2 < y1) ? -1 : 0;

    double tMaxX, tMaxY;
    if (stepX != 0) {
        double nextBoundaryX = (stepX > 0) ? (ix + 1) * cellWidth + minX : ix * cellWidth + minX;
        tMaxX = (nextBoundaryX - x1) / (x2 - x1);
    } else {
        tMaxX = std::numeric_limits<double>::max();
    }

    if (stepY != 0) {
        double nextBoundaryY = (stepY > 0) ? (iy + 1) * cellHeight + minY : iy * cellHeight + minY;
        tMaxY = (nextBoundaryY - y1) / (y2 - y1);
    } else {
        tMaxY = std::numeric_limits<double>::max();
    }

    double tDeltaX = (stepX != 0) ? std::abs(cellWidth / (x2 - x1)) : std::numeric_limits<double>::max();
    double tDeltaY = (stepY != 0) ? std::abs(cellHeight / (y2 - y1)) : std::numeric_limits<double>::max();

    while (true) {
        cells[iy * numCellsX + ix].edgeIndices.push_back(edgeIdx);

        if (ix == iendX && iy == iendY) break;

        if (tMaxX < tMaxY) {
            tMaxX += tDeltaX;
            ix += stepX;
        } else {
            tMaxY += tDeltaY;
            iy += stepY;
        }

        if (ix < 0 || ix >= numCellsX || iy < 0 || iy >= numCellsY) break;
    }
}

void SpatialGrid::removeEdge(int edgeIdx, double x1, double y1, double x2, double y2) {
    if (cells.empty()) return;
    // --- EXACT SAME DDA SETUP AS INSERT EDGE ---
    int ix = static_cast<int>(std::floor((x1 - minX) / cellWidth));
    int iy = static_cast<int>(std::floor((y1 - minY) / cellHeight));
    int iendX = static_cast<int>(std::floor((x2 - minX) / cellWidth));
    int iendY = static_cast<int>(std::floor((y2 - minY) / cellHeight));

    ix = std::clamp(ix, 0, numCellsX - 1);
    iy = std::clamp(iy, 0, numCellsY - 1);
    iendX = std::clamp(iendX, 0, numCellsX - 1);
    iendY = std::clamp(iendY, 0, numCellsY - 1);

    int stepX = (x2 > x1) ? 1 : (x2 < x1) ? -1 : 0;
    int stepY = (y2 > y1) ? 1 : (y2 < y1) ? -1 : 0;

    double tMaxX, tMaxY;
    if (stepX != 0) {
        double nextBoundaryX = (stepX > 0) ? (ix + 1) * cellWidth + minX : ix * cellWidth + minX;
        tMaxX = (nextBoundaryX - x1) / (x2 - x1);
    } else {
        tMaxX = std::numeric_limits<double>::max();
    }

    if (stepY != 0) {
        double nextBoundaryY = (stepY > 0) ? (iy + 1) * cellHeight + minY : iy * cellHeight + minY;
        tMaxY = (nextBoundaryY - y1) / (y2 - y1);
    } else {
        tMaxY = std::numeric_limits<double>::max();
    }

    double tDeltaX = (stepX != 0) ? std::abs(cellWidth / (x2 - x1)) : std::numeric_limits<double>::max();
    double tDeltaY = (stepY != 0) ? std::abs(cellHeight / (y2 - y1)) : std::numeric_limits<double>::max();

    while (true) {
        int cellIndex = iy * numCellsX + ix;
        auto& edgeList = cells[cellIndex].edgeIndices;
        
        auto it = std::find(edgeList.begin(), edgeList.end(), edgeIdx);
        if (it != edgeList.end()) {
            edgeList.erase(it);
        }

        if (ix == iendX && iy == iendY) break;

        if (tMaxX < tMaxY) {
            tMaxX += tDeltaX;
            ix += stepX;
        } else {
            tMaxY += tDeltaY;
            iy += stepY;
        }

        if (ix < 0 || ix >= numCellsX || iy < 0 || iy >= numCellsY) break;
    }
}