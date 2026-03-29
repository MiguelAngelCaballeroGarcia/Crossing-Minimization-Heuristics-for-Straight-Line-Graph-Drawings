#pragma once
#include "logic/Graph.hpp"
#include "logic/PlanarizedGraph.hpp"
#include "logic/SpatialGrid.hpp"
#include <raylib.h>

class GraphVisualizer {
public:
    GraphVisualizer(int screenWidth, int screenHeight, const char* title);
    ~GraphVisualizer();

    void run(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid);

private:
    int width;
    int height;
    
    Camera2D camera;

    // Layer toggles for debugging
    bool showOriginal = false;
    bool showPlanarized = true;
    bool showGrid = true;

    void handleInput();
    void render(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid);
    void drawGridLayer(const SpatialGrid& grid);
};