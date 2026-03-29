#include "GraphVisualizer.hpp"

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

void GraphVisualizer::run(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid) {
    while (!WindowShouldClose()) {
        handleInput();
        render(originalGraph, planarGraph, grid);
    }
}

void GraphVisualizer::handleInput() {
    // --- Layer Toggles ---
    if (IsKeyPressed(KEY_ONE)) showOriginal = !showOriginal;
    if (IsKeyPressed(KEY_TWO)) showPlanarized = !showPlanarized;
    if (IsKeyPressed(KEY_G)) showGrid = !showGrid;

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

void GraphVisualizer::render(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid) {
    BeginDrawing();
    ClearBackground(RAYWHITE);
    BeginMode2D(camera);

    if (showGrid) drawGridLayer(grid);

    if (showOriginal) {
        for (const auto& edge : originalGraph.edges) {
            // Find nodes in the original graph map
            const auto& u = originalGraph.nodes.at(edge.u_id);
            const auto& v = originalGraph.nodes.at(edge.v_id);
            DrawLineEx({(float)u.x, (float)u.y}, {(float)v.x, (float)v.y}, 3.0f, Fade(GRAY, 0.3f));
        }
    }

    if (showPlanarized) {
        // 1. Draw Planar Edges
        for (const auto& [id, pEdge] : planarGraph.edges) {
            const auto& u = planarGraph.nodes.at(pEdge.u_id);
            const auto& v = planarGraph.nodes.at(pEdge.v_id);
            // Draw planar edges in a distinct color (e.g., Dark Blue-ish)
            DrawLineEx({(float)u.x, (float)u.y}, {(float)v.x, (float)v.y}, 2.0f, DARKBLUE);
        }

        // 2. Draw Planar Nodes using DrawCircleV for float precision
        for (const auto& [id, pNode] : planarGraph.nodes) {
            Vector2 pos = {(float)pNode.x, (float)pNode.y};
            if (pNode.type == PlanarizedGraph::NodeType::ORIGINAL) {
                DrawCircleV(pos, 5.0f, BLUE);
            } else { // CROSSING
                DrawCircleV(pos, 6.0f, RED); // Make crossings slightly larger to spot them
            }
        }
    }

    EndMode2D();

    // UI Overlay (Drawn outside of Camera mode so it stays on screen)
    DrawText("Controls:", 10, 10, 20, DARKGRAY);
    DrawText("[1] Toggle Original Graph", 10, 35, 20, showOriginal ? BLACK : LIGHTGRAY);
    DrawText("[2] Toggle Planarized Graph", 10, 60, 20, showPlanarized ? BLACK : LIGHTGRAY);
    DrawText("[G] Toggle Spatial Grid", 10, 85, 20, showGrid ? BLACK : LIGHTGRAY);
    DrawText("Right Click + Drag to Pan | Scroll to Zoom", 10, 115, 20, GRAY);
    
    DrawFPS(GetScreenWidth() - 100, 10);
    EndDrawing();
}