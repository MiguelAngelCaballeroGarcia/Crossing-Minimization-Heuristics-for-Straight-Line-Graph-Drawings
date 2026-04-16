#include "SpatialGrid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

constexpr double kDdaEpsilon = 1e-9;

template <typename Action>
void emitCellIfValid(int col, int row, int numCellsX, int numCellsY, Action&& action) {
    if (col >= 0 && col < numCellsX && row >= 0 && row < numCellsY) {
        action(row * numCellsX + col);
    }
}

bool clipSegmentToGrid(double minX, double maxX, double minY, double maxY,
                       double& x1, double& y1, double& x2, double& y2) {
    const double originalX1 = x1;
    const double originalY1 = y1;
    const double dx = x2 - x1;
    const double dy = y2 - y1;

    double tMin = 0.0;
    double tMax = 1.0;

    auto clip = [&](double p, double q) {
        if (std::abs(p) < kDdaEpsilon) {
            return q >= -kDdaEpsilon;
        }

        const double ratio = q / p;
        if (p < 0.0) {
            if (ratio > tMax) return false;
            if (ratio > tMin) tMin = ratio;
        } else {
            if (ratio < tMin) return false;
            if (ratio < tMax) tMax = ratio;
        }
        return true;
    };

    if (!clip(-dx, x1 - minX)) return false;
    if (!clip(dx, maxX - x1)) return false;
    if (!clip(-dy, y1 - minY)) return false;
    if (!clip(dy, maxY - y1)) return false;

    if (tMax < tMin) return false;

    x1 = originalX1 + (tMin * dx);
    y1 = originalY1 + (tMin * dy);
    x2 = originalX1 + (tMax * dx);
    y2 = originalY1 + (tMax * dy);
    return true;
}

template <typename Action>
void traverseDdaSupercover(const SpatialGrid& grid,
                           double x1, double y1,
                           double x2, double y2,
                           Action&& action) {
    if (grid.getNumCells() == 0) return;

    if (!clipSegmentToGrid(grid.getMinX(), grid.getMaxX(), grid.getMinY(), grid.getMaxY(), x1, y1, x2, y2)) {
        return;
    }

    const double minX = grid.getMinX();
    const double minY = grid.getMinY();
    const double cellWidth = grid.getCellWidth();
    const double cellHeight = grid.getCellHeight();
    const int numCellsX = grid.getNumCellsX();
    const int numCellsY = grid.getNumCellsY();

    const double dx = x2 - x1;
    const double dy = y2 - y1;

    if (std::abs(dx) <= kDdaEpsilon && std::abs(dy) <= kDdaEpsilon) {
        int ix = static_cast<int>(std::floor((x1 - minX) / cellWidth));
        int iy = static_cast<int>(std::floor((y1 - minY) / cellHeight));
        ix = std::clamp(ix, 0, numCellsX - 1);
        iy = std::clamp(iy, 0, numCellsY - 1);
        action(iy * numCellsX + ix);
        return;
    }

    if (std::abs(dx) <= kDdaEpsilon) {
        const double normalizedX = (x1 - minX) / cellWidth;
        const double nearestGridLine = std::round(normalizedX);
        if (std::abs(normalizedX - nearestGridLine) <= kDdaEpsilon) {
            int boundary = static_cast<int>(nearestGridLine);
            int leftCol = boundary - 1;
            int rightCol = boundary;

            if (boundary <= 0) {
                leftCol = -1;
                rightCol = 0;
            } else if (boundary >= numCellsX) {
                leftCol = numCellsX - 1;
                rightCol = numCellsX;
            }

            int rowStart = static_cast<int>(std::floor((std::min(y1, y2) - minY) / cellHeight));
            int rowEnd = static_cast<int>(std::floor((std::max(y1, y2) - minY) / cellHeight));
            rowStart = std::clamp(rowStart, 0, numCellsY - 1);
            rowEnd = std::clamp(rowEnd, 0, numCellsY - 1);

            for (int row = rowStart; row <= rowEnd; ++row) {
                emitCellIfValid(leftCol, row, numCellsX, numCellsY, action);
                emitCellIfValid(rightCol, row, numCellsX, numCellsY, action);
            }
            return;
        }
    }

    if (std::abs(dy) <= kDdaEpsilon) {
        const double normalizedY = (y1 - minY) / cellHeight;
        const double nearestGridLine = std::round(normalizedY);
        if (std::abs(normalizedY - nearestGridLine) <= kDdaEpsilon) {
            int boundary = static_cast<int>(nearestGridLine);
            int topRow = boundary - 1;
            int bottomRow = boundary;

            if (boundary <= 0) {
                topRow = -1;
                bottomRow = 0;
            } else if (boundary >= numCellsY) {
                topRow = numCellsY - 1;
                bottomRow = numCellsY;
            }

            int colStart = static_cast<int>(std::floor((std::min(x1, x2) - minX) / cellWidth));
            int colEnd = static_cast<int>(std::floor((std::max(x1, x2) - minX) / cellWidth));
            colStart = std::clamp(colStart, 0, numCellsX - 1);
            colEnd = std::clamp(colEnd, 0, numCellsX - 1);

            for (int col = colStart; col <= colEnd; ++col) {
                emitCellIfValid(col, topRow, numCellsX, numCellsY, action);
                emitCellIfValid(col, bottomRow, numCellsX, numCellsY, action);
            }
            return;
        }
    }

    int ix = static_cast<int>(std::floor((x1 - minX) / cellWidth));
    int iy = static_cast<int>(std::floor((y1 - minY) / cellHeight));
    int endX = static_cast<int>(std::floor((x2 - minX) / cellWidth));
    int endY = static_cast<int>(std::floor((y2 - minY) / cellHeight));

    ix = std::clamp(ix, 0, numCellsX - 1);
    iy = std::clamp(iy, 0, numCellsY - 1);
    endX = std::clamp(endX, 0, numCellsX - 1);
    endY = std::clamp(endY, 0, numCellsY - 1);

    const int stepX = (dx > 0.0) ? 1 : (dx < 0.0 ? -1 : 0);
    const int stepY = (dy > 0.0) ? 1 : (dy < 0.0 ? -1 : 0);

    const double invDx = (stepX != 0) ? (1.0 / dx) : 0.0;
    const double invDy = (stepY != 0) ? (1.0 / dy) : 0.0;

    double tMaxX = std::numeric_limits<double>::infinity();
    double tMaxY = std::numeric_limits<double>::infinity();
    double tDeltaX = std::numeric_limits<double>::infinity();
    double tDeltaY = std::numeric_limits<double>::infinity();

    if (stepX != 0) {
        const double nextBoundaryX = minX + ((stepX > 0) ? ((ix + 1) * cellWidth) : (ix * cellWidth));
        tMaxX = (nextBoundaryX - x1) * invDx;
        tDeltaX = std::abs(cellWidth * invDx);
    }

    if (stepY != 0) {
        const double nextBoundaryY = minY + ((stepY > 0) ? ((iy + 1) * cellHeight) : (iy * cellHeight));
        tMaxY = (nextBoundaryY - y1) * invDy;
        tDeltaY = std::abs(cellHeight * invDy);
    }

    while (true) {
        emitCellIfValid(ix, iy, numCellsX, numCellsY, action);

        if (ix == endX && iy == endY) {
            break;
        }

        if (std::abs(tMaxX - tMaxY) <= kDdaEpsilon) {
            emitCellIfValid(ix + stepX, iy, numCellsX, numCellsY, action);
            emitCellIfValid(ix, iy + stepY, numCellsX, numCellsY, action);

            tMaxX += tDeltaX;
            tMaxY += tDeltaY;
            ix += stepX;
            iy += stepY;
        } else if (tMaxX < tMaxY) {
            tMaxX += tDeltaX;
            ix += stepX;
        } else {
            tMaxY += tDeltaY;
            iy += stepY;
        }

        if (ix < 0 || ix >= numCellsX || iy < 0 || iy >= numCellsY) {
            break;
        }
    }
}

} // namespace

SpatialGrid::SpatialGrid()
    : minX(0), maxX(0), minY(0), maxY(0),
      cellWidth(1), cellHeight(1), numCellsX(0), numCellsY(0) {}

SpatialGrid::SpatialGrid(double minX, double maxX, double minY, double maxY, int expectedTotalNodes)
    : minX(minX), maxX(maxX), minY(minY), maxY(maxY)
{
    double width = maxX - minX;
    double height = maxY - minY;

    if (width == 0.0) width = 1.0;
    if (height == 0.0) height = 1.0;

    double targetTotalCells = static_cast<double>(expectedTotalNodes) / 5.0;
    if (targetTotalCells < 1.0) targetTotalCells = 1.0;

    const double aspectRatio = width / height;
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

    ix = std::clamp(ix, 0, numCellsX - 1);
    iy = std::clamp(iy, 0, numCellsY - 1);

    return iy * numCellsX + ix;
}

void SpatialGrid::insertNode(int nodeIdx, double x, double y) {
    if (cells.empty()) return;
    cells[getCellIndex(x, y)].nodeIndices.push_back(nodeIdx);
}

void SpatialGrid::removeNode(int nodeIdx, double x, double y) {
    if (cells.empty()) return;

    auto& nodeList = cells[getCellIndex(x, y)].nodeIndices;
    auto it = std::find(nodeList.begin(), nodeList.end(), nodeIdx);
    if (it != nodeList.end()) {
        std::swap(*it, nodeList.back());
        nodeList.pop_back();
    }
}

void SpatialGrid::insertEdge(int edgeIdx, double x1, double y1, double x2, double y2) {
    traverseDdaSupercover(*this, x1, y1, x2, y2, [&](int cellIndex) {
        auto& edgeList = cells[cellIndex].edgeIndices;
        if (std::find(edgeList.begin(), edgeList.end(), edgeIdx) == edgeList.end()) {
            edgeList.push_back(edgeIdx);
        }
    });
}

void SpatialGrid::removeEdge(int edgeIdx, double x1, double y1, double x2, double y2) {
    traverseDdaSupercover(*this, x1, y1, x2, y2, [&](int cellIndex) {
        auto& edgeList = cells[cellIndex].edgeIndices;
        auto it = std::find(edgeList.begin(), edgeList.end(), edgeIdx);
        while (it != edgeList.end()) {
            std::swap(*it, edgeList.back());
            edgeList.pop_back();
            it = std::find(edgeList.begin(), edgeList.end(), edgeIdx);
        }
    });
}