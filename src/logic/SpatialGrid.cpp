#include "SpatialGrid.hpp"

#include <cmath>
#include <algorithm>
#include <limits>


int SpatialGrid::getCellIndex(double x, double y) const {
    // 1. Calculate relative offset
    int ix = static_cast<int>(std::floor((x - minX) / cellWidth));
    int iy = static_cast<int>(std::floor((y - minY) / cellHeight));

    // 2. Clamp to ensure we stay inside the vector bounds
    // (Prevents crashes if a node is exactly on the maxX/maxY boundary)
    if (ix < 0) ix = 0;
    if (ix >= numCellsX) ix = numCellsX - 1;
    if (iy < 0) iy = 0;
    if (iy >= numCellsY) iy = numCellsY - 1;

    return iy * numCellsX + ix;
}

// NEW Bulk-Loading Constructor
SpatialGrid::SpatialGrid(const Graph& graph) {
    // 1. Find the exact bounding box of the graph
    minX = std::numeric_limits<double>::max();
    minY = std::numeric_limits<double>::max();
    maxX = std::numeric_limits<double>::lowest();
    maxY = std::numeric_limits<double>::lowest();

    for (const auto& node : graph.nodes) {
        if (node.x < minX) minX = node.x;
        if (node.x > maxX) maxX = node.x;
        if (node.y < minY) minY = node.y;
        if (node.y > maxY) maxY = node.y;
    }

    // 2. Setup Dimensions (Using your excellent 5-node logic)
    int totalNodes = graph.nodes.size();
    double width = maxX - minX;
    double height = maxY - minY;
    
    double targetTotalCells = static_cast<double>(totalNodes) / 5.0;
    if (targetTotalCells < 1) targetTotalCells = 1;

    double aspectRatio = width / height;
    numCellsX = static_cast<int>(std::sqrt(targetTotalCells * aspectRatio));
    numCellsX = std::max(1, numCellsX);
    numCellsY = static_cast<int>(std::ceil(targetTotalCells / numCellsX));
    numCellsY = std::max(1, numCellsY);

    cellWidth = width / numCellsX;
    cellHeight = height / numCellsY;

    cells.resize(numCellsX * numCellsY);

    // 3. Bulk Insert Nodes
    for (const auto& node : graph.nodes) {
        insertNode(node.id, node.x, node.y);
    }

    // 4. Bulk Insert Edges
    for (const auto& edge : graph.edges) {
        const auto& u = graph.nodes[edge.u_id]; // Assuming edge has u, v indices
        const auto& v = graph.nodes[edge.v_id];
        insertEdge(edge.id, u.x, u.y, v.x, v.y);
    }
}

void SpatialGrid::insertNode(int nodeIdx, double x, double y) {
    // Determine which cell the node belongs to
    int cellIndex = getCellIndex(x, y);
    // Insert the node into the appropriate cell
    cells[cellIndex].nodeIndices.push_back(nodeIdx);
}

void SpatialGrid::insertEdge(int edgeIdx, double x1, double y1, double x2, double y2) {
    // 1. Identify the starting cell coordinates
    int ix = static_cast<int>(std::floor((x1 - minX) / cellWidth));
    int iy = static_cast<int>(std::floor((y1 - minY) / cellHeight));
    
    // 2. Identify the ending cell coordinates
    int iendX = static_cast<int>(std::floor((x2 - minX) / cellWidth));
    int iendY = static_cast<int>(std::floor((y2 - minY) / cellHeight));

    // Clamp coordinates to ensure we stay within the grid bounds
    ix = std::clamp(ix, 0, numCellsX - 1);
    iy = std::clamp(iy, 0, numCellsY - 1);
    iendX = std::clamp(iendX, 0, numCellsX - 1);
    iendY = std::clamp(iendY, 0, numCellsY - 1);

    // 3. Determine the direction of the step (+1 or -1)
    int stepX = (x2 > x1) ? 1 : (x2 < x1) ? -1 : 0;
    int stepY = (y2 > y1) ? 1 : (y2 < y1) ? -1 : 0;

    // 4. Calculate the distance to the next horizontal/vertical cell boundary
    // tMax is the value of 't' (from 0 to 1 along the segment) at the next boundary
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

    // 5. Calculate how far 't' changes when we cross exactly one cell width/height
    double tDeltaX = (stepX != 0) ? std::abs(cellWidth / (x2 - x1)) : std::numeric_limits<double>::max();
    double tDeltaY = (stepY != 0) ? std::abs(cellHeight / (y2 - y1)) : std::numeric_limits<double>::max();

    // 6. Traversal Loop
    while (true) {
        // Register the edge in the current cell
        cells[iy * numCellsX + ix].edgeIndices.push_back(edgeIdx);

        // If we've reached the target cell, we are done
        if (ix == iendX && iy == iendY) break;

        // Move to the next closest boundary (X or Y)
        if (tMaxX < tMaxY) {
            tMaxX += tDeltaX;
            ix += stepX;
        } else {
            tMaxY += tDeltaY;
            iy += stepY;
        }

        // Safety check to avoid infinite loops or out-of-bounds
        if (ix < 0 || ix >= numCellsX || iy < 0 || iy >= numCellsY) break;
    }
}

// NEW: Remove Edge logic for incremental updates
void SpatialGrid::removeEdge(int edgeIdx, double x1, double y1, double x2, double y2) {
    // This uses the EXACT same DDA traversal logic as insertEdge
    // I will abbreviate the DDA setup for space, but you copy the ix, iy, stepX, etc. from insertEdge

    // 1. Identify the starting cell coordinates
    int ix = static_cast<int>(std::floor((x1 - minX) / cellWidth));
    int iy = static_cast<int>(std::floor((y1 - minY) / cellHeight));
    
    // 2. Identify the ending cell coordinates
    int iendX = static_cast<int>(std::floor((x2 - minX) / cellWidth));
    int iendY = static_cast<int>(std::floor((y2 - minY) / cellHeight));

    // Clamp coordinates to ensure we stay within the grid bounds
    ix = std::clamp(ix, 0, numCellsX - 1);
    iy = std::clamp(iy, 0, numCellsY - 1);
    iendX = std::clamp(iendX, 0, numCellsX - 1);
    iendY = std::clamp(iendY, 0, numCellsY - 1);

    // 3. Determine the direction of the step (+1 or -1)
    int stepX = (x2 > x1) ? 1 : (x2 < x1) ? -1 : 0;
    int stepY = (y2 > y1) ? 1 : (y2 < y1) ? -1 : 0;

    // 4. Calculate the distance to the next horizontal/vertical cell boundary
    // tMax is the value of 't' (from 0 to 1 along the segment) at the next boundary
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

    // 5. Calculate how far 't' changes when we cross exactly one cell width/height
    double tDeltaX = (stepX != 0) ? std::abs(cellWidth / (x2 - x1)) : std::numeric_limits<double>::max();
    double tDeltaY = (stepY != 0) ? std::abs(cellHeight / (y2 - y1)) : std::numeric_limits<double>::max();

    while (true) {
        int cellIndex = iy * numCellsX + ix;
        auto& edgeList = cells[cellIndex].edgeIndices;
        
        // Find and erase the edge ID from this cell
        auto it = std::find(edgeList.begin(), edgeList.end(), edgeIdx);
        if (it != edgeList.end()) {
            edgeList.erase(it);
        }

        if (ix == iendX && iy == iendY) break;

        // Move to the next closest boundary (X or Y)
        if (tMaxX < tMaxY) {
            tMaxX += tDeltaX;
            ix += stepX;
        } else {
            tMaxY += tDeltaY;
            iy += stepY;
        }

        // Safety check to avoid infinite loops or out-of-bounds
        if (ix < 0 || ix >= numCellsX || iy < 0 || iy >= numCellsY) break;
    }
}