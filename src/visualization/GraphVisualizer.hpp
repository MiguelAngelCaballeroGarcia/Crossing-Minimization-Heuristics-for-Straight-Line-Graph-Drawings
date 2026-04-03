#pragma once
#include "logic/Graph.hpp"
#include "logic/PlanarizedGraph.hpp"
#include "logic/SpatialGrid.hpp"
#include "../heuristics/RelocationManager.hpp" // Added
#include <raylib.h>
#include <vector>

class GraphVisualizer {
public:
    GraphVisualizer(int screenWidth, int screenHeight, const char* title);
    ~GraphVisualizer();

    void run(const Graph& originalGraph, PlanarizedGraph& planarGraph, const SpatialGrid& grid);

private:
    int width;
    int height;
    Camera2D camera;

    // Layer toggles
    bool showOriginal = false;
    bool showPlanarized = true;
    bool showGrid = true;
    bool showROI = true; // New toggle

    // Relocation Debug Data
    int activeNodeId = -1;
    RegionOfInterest currentROI;
    std::vector<LocalSegment> currentLocalSegments;

    void handleInput(PlanarizedGraph& planarGraph);
    void render(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid);
    void drawGridLayer(const SpatialGrid& grid);
    void drawROILayer(); // New helper
    int countCrossingNodes(const PlanarizedGraph& planarGraph) const;
};