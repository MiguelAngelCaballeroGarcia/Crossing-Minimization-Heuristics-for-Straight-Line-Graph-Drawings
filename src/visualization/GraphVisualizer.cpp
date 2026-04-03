#include "GraphVisualizer.hpp"
#include <algorithm>

GraphVisualizer::GraphVisualizer(int screenWidth, int screenHeight, const char* title)
    : width(screenWidth), height(screenHeight) {
    InitWindow(width, height, title);
    SetTargetFPS(60);

    // Initialize 2D Camera for panning and zooming
    camera = { 0 };
    camera.target = { (float)width / 2, (float)height / 2 };
    camera.offset = { (float)width / 2, (float)height / 2 };
    camera.rotation = 0.0f;
    camera.zoom = 1.0f;
}

GraphVisualizer::~GraphVisualizer() {
    CloseWindow();
}

void GraphVisualizer::run(const Graph& originalGraph, PlanarizedGraph& planarGraph, const SpatialGrid& grid) {
    while (!WindowShouldClose()) {
        handleInput(planarGraph);
        render(originalGraph, planarGraph, grid);
    }
}

void GraphVisualizer::handleInput(PlanarizedGraph& planarGraph) {
    // --- Layer Toggles ---
    if (IsKeyPressed(KEY_ONE)) showOriginal = !showOriginal;
    if (IsKeyPressed(KEY_TWO)) showPlanarized = !showPlanarized;
    if (IsKeyPressed(KEY_G)) showGrid = !showGrid;
    if (IsKeyPressed(KEY_H)) showROI = !showROI; // H for "Heuristic/Highlight"

    // --- TRIGGER RELOCATION STEP (PHASE 1) ---
    if (IsKeyPressed(KEY_R)) {
        RelocationManager manager(planarGraph);
        
        activeNodeId = manager.selectVariableNode();
        if (activeNodeId != -1) {
            currentROI = manager.calculateROI(activeNodeId);
            currentLocalSegments = manager.extractAndClipGeometry(currentROI, activeNodeId);
        }
    }

    // --- Camera Panning (Right Mouse Button) ---
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
        Vector2 delta = GetMouseDelta();
        delta.x *= -1.0f / camera.zoom;
        delta.y *= -1.0f / camera.zoom;
        camera.target.x += delta.x;
        camera.target.y += delta.y;
    }

    // --- Camera Zooming (Mouse Wheel) ---
    float wheel = GetMouseWheelMove();
    if (wheel != 0) {
        // Get the world point that is under the mouse
        Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), camera);
        
        // Set the offset to where the mouse is
        camera.offset = GetMousePosition();
        
        // Set the target to match, so that the camera zooms around the mouse position
        camera.target = mouseWorldPos;
        
        // Apply zoom
        const float zoomIncrement = 0.125f;
        camera.zoom += (wheel * zoomIncrement);
        if (camera.zoom < 0.1f) camera.zoom = 0.1f;
    }
}

void GraphVisualizer::drawROILayer() {
    if (activeNodeId == -1 || !showROI) return;

    // 1. Draw the ROI Boundary Box (The "Frozen" area)
    Rectangle roiRect = {
        (float)currentROI.minX, 
        (float)currentROI.minY, 
        (float)(currentROI.maxX - currentROI.minX), 
        (float)(currentROI.maxY - currentROI.minY)
    };
    DrawRectangleLinesEx(roiRect, 4.0f, YELLOW);
    DrawRectangleRec(roiRect, Fade(YELLOW, 0.1f));

    // 2. Draw the Clipped Local Segments
    for (const auto& seg : currentLocalSegments) {
        Color segmentColor = (seg.originalEdgeId == -1) ? RED : ORANGE;
        float thickness = (seg.originalEdgeId == -1) ? 3.0f : 2.0f;
        
        DrawLineEx(
            {(float)seg.x1, (float)seg.y1}, 
            {(float)seg.x2, (float)seg.y2}, 
            thickness, 
            segmentColor
        );
    }
}

void GraphVisualizer::drawGridLayer(const SpatialGrid& grid) {
    // Use float for everything to avoid snapping to integers
    float minX = (float)grid.getMinX();
    float minY = (float)grid.getMinY();
    float maxX = (float)grid.getMaxX();
    float maxY = (float)grid.getMaxY();
    float cWidth = (float)grid.getCellWidth();
    float cHeight = (float)grid.getCellHeight();

    // Vertical Lines
    for (int i = 0; i <= grid.getNumCellsX(); i++) {
        float x = minX + (i * cWidth);
        DrawLineEx({x, minY}, {x, maxY}, 1.0f, Fade(LIGHTGRAY, 0.5f));
    }

    // Horizontal Lines
    for (int i = 0; i <= grid.getNumCellsY(); i++) {
        float y = minY + (i * cHeight);
        DrawLineEx({minX, y}, {maxX, y}, 1.0f, Fade(LIGHTGRAY, 0.5f));
    }
    
    DrawRectangleLinesEx({minX, minY, maxX - minX, maxY - minY}, 2.0f, GREEN);
}

int GraphVisualizer::countCrossingNodes(const PlanarizedGraph& planarGraph) const {
    return static_cast<int>(std::count_if(
        planarGraph.nodes.begin(),
        planarGraph.nodes.end(),
        [](const auto& pair) {
            return pair.second.type == PlanarizedGraph::NodeType::CROSSING;
        }
    ));
}

void GraphVisualizer::render(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid) {
    BeginDrawing();
    ClearBackground(RAYWHITE);
    BeginMode2D(camera);

    if (showGrid) drawGridLayer(grid);
    
    // Draw ROI and Local Segments behind the nodes but on top of grid
    drawROILayer();

    if (showOriginal) {
        for (const auto& edge : originalGraph.edges) {
            const auto& u = originalGraph.nodes.at(edge.u_id);
            const auto& v = originalGraph.nodes.at(edge.v_id);
            DrawLineEx({(float)u.x, (float)u.y}, {(float)v.x, (float)v.y}, 3.0f, Fade(GRAY, 0.2f));
        }
    }

    if (showPlanarized) {
        // Draw Edges
        for (const auto& [id, pEdge] : planarGraph.edges) {
            const auto& u = planarGraph.nodes.at(pEdge.u_id);
            const auto& v = planarGraph.nodes.at(pEdge.v_id);
            DrawLineEx({(float)u.x, (float)u.y}, {(float)v.x, (float)v.y}, 1.0f, DARKBLUE);
        }

        // Draw Nodes
        for (const auto& [id, pNode] : planarGraph.nodes) {
            Vector2 pos = {(float)pNode.x, (float)pNode.y};
            
            // Highlight the active node being relocated
            if (id == activeNodeId) {
                DrawCircleV(pos, 6.0f, GOLD);
                DrawCircleLines(pos.x, pos.y, 10.0f, BLACK);
            } else if (pNode.type == PlanarizedGraph::NodeType::ORIGINAL) {
                DrawCircleV(pos, 4.0f, BLUE);
            } else {
                DrawCircleV(pos, 2.0f, RED);
            }
        }
    }

    EndMode2D();

    // UI Overlay (Drawn outside of Camera mode so it stays on screen)
    DrawText("Controls:", 10, 10, 20, DARKGRAY);
    DrawText("[R] Randomly Select Node & Extract ROI", 10, 35, 20, MAROON);
    DrawText("[H] Toggle ROI Visualization", 10, 60, 20, showROI ? BLACK : LIGHTGRAY);
    DrawText("[1] Toggle Original Graph", 10, 85, 20, showOriginal ? BLACK : LIGHTGRAY);
    DrawText("[2] Toggle Planarized Graph", 10, 110, 20, showPlanarized ? BLACK : LIGHTGRAY);
    DrawText("[G] Toggle Spatial Grid", 10, 135, 20, showGrid ? BLACK : LIGHTGRAY);
    DrawText(TextFormat("Crossings: %d", countCrossingNodes(planarGraph)), 10, 160, 20, MAROON);
    DrawText("Right Click + Drag to Pan | Scroll to Zoom", 10, 185, 20, GRAY);
    
    DrawFPS(GetScreenWidth() - 100, 10);
    EndDrawing();
}