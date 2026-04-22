#pragma once
#include "logic/Graph.hpp"
#include "logic/PlanarizedGraph.hpp"
#include "logic/SpatialGrid.hpp"
#include "../heuristics/RelocationManager.hpp" // Added
#include <raylib.h>
#include <vector>
#include <string>

class GraphVisualizer {
public:
    GraphVisualizer(int screenWidth, int screenHeight, const char* title);
    ~GraphVisualizer();

    void run(const Graph& originalGraph, PlanarizedGraph& planarGraph, const SpatialGrid& grid);

private:
    int width;
    int height;
    Camera2D camera;
    
    const Graph* originalGraphPtr = nullptr;  // Store original graph for exact crossing detection

    // Layer toggles
    bool showOriginal = false;
    bool showPlanarized = true;
    bool showGrid = true;
    bool showROI = true; // New toggle
    bool showNodeIds = false;
    bool showScaleRuler = false;
    ROIBoundaryMethod roiBoundaryMethod = ROIBoundaryMethod::NEIGHBORS_INSIDE;

    // Relocation Debug Data
    int activeNodeId = -1;
    RegionOfInterest currentROI;
    std::vector<LocalSegment> currentLocalSegments;
    LocalRegionAnalysis currentAnalysis;

    bool iterationBoxActive = false;
    std::string iterationInput = "100";
    int lastBatchIterations = 0;
    int lastBatchMoves = 0;
    bool lastBatchExplorationUsed = false;
    int lastBatchExactCrossings = -1;
    int lastBatchNodeBasedCrossings = -1;
    int lastBatchCrossingDiff = 0;

    void handleInput(PlanarizedGraph& planarGraph);
    void render(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid);
    void drawGridLayer(const SpatialGrid& grid);
    void drawScaleRuler(const SpatialGrid& grid);
    void drawROILayer(); // New helper
    void drawDualGraphPanel();
    int countCrossingNodes(const PlanarizedGraph& planarGraph) const;
};