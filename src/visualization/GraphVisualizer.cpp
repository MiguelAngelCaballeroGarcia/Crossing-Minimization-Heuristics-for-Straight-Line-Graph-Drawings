#include "GraphVisualizer.hpp"
#include "geometry/IntersectionDetector.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

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
    originalGraphPtr = &originalGraph;
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
    if (IsKeyPressed(KEY_N)) showNodeIds = !showNodeIds;
    if (IsKeyPressed(KEY_K)) showScaleRuler = !showScaleRuler;

    const Rectangle iterationBox = {(float)width - 250.0f, 10.0f, 180.0f, 32.0f};
    const Rectangle roiBoundaryButton = {(float)width - 250.0f, 170.0f, 220.0f, 32.0f};
    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
        const Vector2 mousePos = GetMousePosition();

        if (CheckCollisionPointRec(mousePos, roiBoundaryButton)) {
            switch (roiBoundaryMethod) {
                case ROIBoundaryMethod::NEIGHBORS_INSIDE:
                    roiBoundaryMethod = ROIBoundaryMethod::FULL_GRID;
                    break;
                case ROIBoundaryMethod::FULL_GRID:
                    roiBoundaryMethod = ROIBoundaryMethod::LOCAL_3X3_AROUND_NODE_CELL;
                    break;
                case ROIBoundaryMethod::LOCAL_3X3_AROUND_NODE_CELL:
                default:
                    roiBoundaryMethod = ROIBoundaryMethod::NEIGHBORS_INSIDE;
                    break;
            }

            iterationBoxActive = false;

            if (activeNodeId != -1) {
                RelocationManager manager(planarGraph);
                manager.setROIBoundaryMethod(roiBoundaryMethod);
                currentROI = manager.calculateROI(activeNodeId);
                currentAnalysis = manager.analyzeLocalRegions(currentROI, activeNodeId);
                currentLocalSegments = currentAnalysis.localGeometry;
            }
        } else {
            iterationBoxActive = CheckCollisionPointRec(mousePos, iterationBox);
        }
    }

    if (iterationBoxActive) {
        int key = GetCharPressed();
        while (key > 0) {
            if (key >= '0' && key <= '9' && iterationInput.size() < 7) {
                iterationInput.push_back(static_cast<char>(key));
            }
            key = GetCharPressed();
        }

        if (IsKeyPressed(KEY_BACKSPACE) && !iterationInput.empty()) {
            iterationInput.pop_back();
        }

        if (IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_KP_ENTER)) {
            int iterations = 0;
            if (!iterationInput.empty()) {
                iterations = std::atoi(iterationInput.c_str());
            }
            iterations = std::clamp(iterations, 0, 100000);

            lastBatchIterations = iterations;
            lastBatchMoves = 0;
            lastBatchExplorationUsed = false;
            lastBatchExactCrossings = -1;
            lastBatchNodeBasedCrossings = -1;
            lastBatchCrossingDiff = 0;

            if (iterations > 0) {
                RelocationManager manager(planarGraph);
                manager.setROIBoundaryMethod(roiBoundaryMethod);

                for (int i = 0; i < iterations; ++i) {
                    RelocationStepResult step = manager.performRelocationStep();
                    if (step.moved) {
                        lastBatchMoves++;
                    }
                    if (step.usedExplorationMove) {
                        lastBatchExplorationUsed = true;
                    }
                }

                if (originalGraphPtr) {
                    Graph currentGraph = *originalGraphPtr;

                    for (size_t i = 0; i < currentGraph.nodes.size(); ++i) {
                        const int planarNodeId = currentGraph.nodes[i].id;
                        if (!planarGraph.hasNode(planarNodeId)) continue;
                        const auto& pNode = planarGraph.getNode(planarNodeId);
                        currentGraph.nodes[i].x = pNode.x;
                        currentGraph.nodes[i].y = pNode.y;
                    }

                    auto [minX, minY, maxX, maxY] = currentGraph.getBounds();
                    SpatialGrid tempGrid(minX, maxX, minY, maxY, currentGraph.nodes.size());

                    for (const auto& edge : currentGraph.edges) {
                        if (edge.u_id < 0 || edge.v_id < 0 ||
                            edge.u_id >= static_cast<int>(currentGraph.nodes.size()) ||
                            edge.v_id >= static_cast<int>(currentGraph.nodes.size())) {
                            continue;
                        }

                        const auto& u = currentGraph.nodes[edge.u_id];
                        const auto& v = currentGraph.nodes[edge.v_id];
                        tempGrid.insertEdge(edge.id, u.x, u.y, v.x, v.y);
                    }

                    const auto exactIntersections = IntersectionDetector::findIntersections(currentGraph, tempGrid);
                    lastBatchExactCrossings = static_cast<int>(exactIntersections.size());
                    lastBatchNodeBasedCrossings = planarGraph.countTotalCrossings();
                    lastBatchCrossingDiff = lastBatchExactCrossings - lastBatchNodeBasedCrossings;

                    std::cerr << "[Batch Complete] Iterations: " << iterations
                              << " | Moves: " << lastBatchMoves
                              << " | Exact crossings: " << lastBatchExactCrossings
                              << " | Node-based crossings: " << lastBatchNodeBasedCrossings
                              << " | Diff: " << lastBatchCrossingDiff << "\n";
                }
            }

            // Batch mode renders only final state, so clear per-step ROI overlays.
            activeNodeId = -1;
            currentLocalSegments.clear();
            currentAnalysis = LocalRegionAnalysis{};
        }
    }

    // --- TRIGGER RELOCATION STEP (PHASE 1) ---
    if (IsKeyPressed(KEY_R)) {
        RelocationManager manager(planarGraph);
        manager.setROIBoundaryMethod(roiBoundaryMethod);
        
        activeNodeId = manager.selectVariableNode();
        if (activeNodeId != -1) {
            currentROI = manager.calculateROI(activeNodeId);
            currentAnalysis = manager.analyzeLocalRegions(currentROI, activeNodeId);
            currentLocalSegments = currentAnalysis.localGeometry;
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

    // 2. Draw the Clipped Local Segments with type-based colors
    for (const auto& seg : currentLocalSegments) {
        Color segmentColor = ORANGE;
        float thickness = 2.0f;

        if (seg.type == LocalSegmentType::RAY) {
            segmentColor = RED;
            thickness = 3.0f;
        } else if (seg.type == LocalSegmentType::ROI_BORDER) {
            segmentColor = {139, 69, 19, 255}; // Brown
            thickness = 3.0f;
        }
        // ORIGINAL_SUBSEGMENT uses ORANGE
        
        DrawLineEx(
            {(float)seg.x1, (float)seg.y1}, 
            {(float)seg.x2, (float)seg.y2}, 
            thickness, 
            segmentColor
        );
    }

    // 3. Draw face weights (excluding outer face)
    for (const auto& face : currentAnalysis.dualGraph.faces) {
        if (face.isOuter || face.vertices.empty()) continue;
        if (face.id >= static_cast<int>(currentAnalysis.faceWeights.size())) continue;

        const double w = currentAnalysis.faceWeights[face.id];
        if (!std::isfinite(w)) continue;

        double cx = 0.0;
        double cy = 0.0;
        for (const auto& v : face.vertices) {
            cx += v.x;
            cy += v.y;
        }
        cx /= static_cast<double>(face.vertices.size());
        cy /= static_cast<double>(face.vertices.size());

        // Color based on whether it's the source face
        Color faceColor = (face.id == currentAnalysis.sourceFaceId) ? BLUE : DARKGREEN;
        
        // Draw weight label
        DrawText(TextFormat("F%d: %.0f", face.id, w), 
            static_cast<int>(cx), static_cast<int>(cy), 9, faceColor);
    }
}

void GraphVisualizer::drawDualGraphPanel() {
    if (activeNodeId == -1 || !showROI) return;
    if (currentAnalysis.dualGraph.faces.empty()) return;

    const int panelX = 10;
    const int panelY = 240;
    const int panelW = std::min(640, GetScreenWidth() - 20);
    const int panelH = std::min(460, GetScreenHeight() - panelY - 10);

    DrawRectangle(panelX, panelY, panelW, panelH, Fade(LIGHTGRAY, 0.5f));
    DrawRectangleLines(panelX, panelY, panelW, panelH, DARKGRAY);
    DrawText("Dual Graph (ROI)  label: F:w<local> d<global>", panelX + 8, panelY + 8, 18, BLACK);

    const double roiW = std::max(1e-9, currentROI.maxX - currentROI.minX);
    const double roiH = std::max(1e-9, currentROI.maxY - currentROI.minY);
    const int margin = 30;

    auto toPanel = [&](double x, double y) {
        const double nx = (x - currentROI.minX) / roiW;
        const double ny = (y - currentROI.minY) / roiH;
        const float px = static_cast<float>(panelX + margin + nx * (panelW - 2 * margin));
        const float py = static_cast<float>(panelY + margin + ny * (panelH - 2 * margin));
        return Vector2{px, py};
    };

    std::vector<Vector2> centers(currentAnalysis.dualGraph.faces.size(), {0.0f, 0.0f});
    for (const auto& face : currentAnalysis.dualGraph.faces) {
        if (face.vertices.empty()) continue;
        double cx = 0.0;
        double cy = 0.0;
        for (const auto& v : face.vertices) {
            cx += v.x;
            cy += v.y;
        }
        cx /= static_cast<double>(face.vertices.size());
        cy /= static_cast<double>(face.vertices.size());
        centers[face.id] = toPanel(cx, cy);
    }

    for (const auto& e : currentAnalysis.dualTreeEdges) {
        if (e.faceA < 0 || e.faceB < 0) continue;
        if (e.faceA >= static_cast<int>(currentAnalysis.dualGraph.faces.size()) ||
            e.faceB >= static_cast<int>(currentAnalysis.dualGraph.faces.size())) continue;
        const auto& aFace = currentAnalysis.dualGraph.faces[e.faceA];
        const auto& bFace = currentAnalysis.dualGraph.faces[e.faceB];
        if (aFace.isOuter || bFace.isOuter) continue;

        Color c = (e.boundaryType == LocalSegmentType::RAY) ? MAGENTA : DARKBLUE;
        DrawLineEx(centers[e.faceA], centers[e.faceB], 3.0f, c);
    }

    for (const auto& face : currentAnalysis.dualGraph.faces) {
        if (face.isOuter) continue;
        const auto p = centers[face.id];
        const Color nodeColor = (face.id == currentAnalysis.sourceFaceId) ? BLUE : BLACK;
        DrawCircleV(p, 6.0f, nodeColor);

        if (face.id >= 0 && face.id < static_cast<int>(currentAnalysis.faceWeights.size()) &&
            std::isfinite(currentAnalysis.faceWeights[face.id])) {
            const int globalDelta = (face.id < static_cast<int>(currentAnalysis.faceGlobalCrossingDelta.size()))
                ? currentAnalysis.faceGlobalCrossingDelta[face.id]
                : 0;
            DrawText(TextFormat("%d:w%.0f d%d", face.id, currentAnalysis.faceWeights[face.id], globalDelta),
                     static_cast<int>(p.x) + 6,
                     static_cast<int>(p.y) - 8,
                     16,
                     nodeColor);
        }
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

void GraphVisualizer::drawScaleRuler(const SpatialGrid& grid) {
    if (!showScaleRuler) return;

    const float minX = static_cast<float>(grid.getMinX());
    const float minY = static_cast<float>(grid.getMinY());
    const float maxX = static_cast<float>(grid.getMaxX());
    const float maxY = static_cast<float>(grid.getMaxY());

    const float graphW = std::max(1.0f, maxX - minX);
    const float graphH = std::max(1.0f, maxY - minY);
    const float marginX = 0.03f * graphW;
    const float marginY = 0.06f * graphH;

    const float originX = minX + marginX;
    const float originY = minY + marginY;

    // World-space reference: ruler length is exactly 100 graph units.
    const float rulerLength = 100.0f;
    DrawLineEx({originX, originY}, {originX + rulerLength, originY}, 2.5f, BLACK);

    for (int i = 0; i <= 100; ++i) {
        const float x = originX + static_cast<float>(i);
        float tickLen = 3.0f;
        Color tickColor = DARKGRAY;

        if (i % 10 == 0) {
            tickLen = 10.0f;
            tickColor = MAROON;
        } else if (i % 5 == 0) {
            tickLen = 6.0f;
            tickColor = GRAY;
        }

        DrawLineEx({x, originY - tickLen}, {x, originY + tickLen}, 1.0f, tickColor);

        if (i % 10 == 0) {
            DrawText(TextFormat("%d", i), static_cast<int>(x - 6.0f), static_cast<int>(originY + 14.0f), 10, MAROON);
        }
    }

    DrawText("Scale ruler in graph units (0..100)",
             static_cast<int>(originX),
             static_cast<int>(originY - 24.0f),
             12,
             BLACK);
}

int GraphVisualizer::countCrossingNodes(const PlanarizedGraph& planarGraph) const {
    int count = 0;
    planarGraph.forEachNode([&count](int, const PlanarizedGraph::PlanarNode& node) {
        if (node.type == PlanarizedGraph::NodeType::CROSSING) {
            ++count;
        }
    });
    return count;
}

void GraphVisualizer::render(const Graph& originalGraph, const PlanarizedGraph& planarGraph, const SpatialGrid& grid) {
    BeginDrawing();
    ClearBackground(RAYWHITE);
    BeginMode2D(camera);

    if (showGrid) drawGridLayer(grid);
    drawScaleRuler(grid);
    
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
        planarGraph.forEachEdge([&](int id, const PlanarizedGraph::PlanarEdge& pEdge) { // 'id' is the edge ID
            const auto& u = planarGraph.getNode(pEdge.u_id);
            const auto& v = planarGraph.getNode(pEdge.v_id);
            
            Vector2 start = {(float)u.x, (float)u.y};
            Vector2 end = {(float)v.x, (float)v.y};

            // Draw the edge line
            DrawLineEx(start, end, 1.0f, DARKBLUE);

            // NEW: Draw the Edge ID if toggled
            if (showNodeIds) {
                // Calculate the midpoint of the edge
                int midX = (int)((start.x + end.x) / 2.0f);
                int midY = (int)((start.y + end.y) / 2.0f);

                // Draw the ID (using a distinct color like DARKGREEN or GRAY)
                DrawText(TextFormat("E%d", pEdge.original_edge_id), midX, midY, 10, DARKGREEN);
            }
        });

        // Draw Nodes
        planarGraph.forEachNode([&](int id, const PlanarizedGraph::PlanarNode& pNode) {
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

            if (showNodeIds) {
                const Color labelColor = (pNode.type == PlanarizedGraph::NodeType::ORIGINAL) ? DARKBLUE : MAROON;
                DrawText(TextFormat("%d", id),
                         static_cast<int>(pos.x) + 6,
                         static_cast<int>(pos.y) - 14,
                         12,
                         labelColor);
            }
        });
    }

    EndMode2D();

    drawDualGraphPanel();

    // UI Overlay (Drawn outside of Camera mode so it stays on screen)
    DrawText("Controls:", 10, 10, 20, DARKGRAY);
    DrawText("[R] Randomly Select Node & Analyze ROI", 10, 35, 20, MAROON);
    DrawText("[H] Toggle ROI Visualization", 10, 60, 20, showROI ? BLACK : LIGHTGRAY);
    DrawText("[1] Toggle Original Graph", 10, 85, 20, showOriginal ? BLACK : LIGHTGRAY);
    DrawText("[2] Toggle Planarized Graph", 10, 110, 20, showPlanarized ? BLACK : LIGHTGRAY);
    DrawText("[G] Toggle Spatial Grid", 10, 135, 20, showGrid ? BLACK : LIGHTGRAY);
    DrawText("[N] Toggle Node IDs", 10, 160, 20, showNodeIds ? BLACK : LIGHTGRAY);
    DrawText("[K] Toggle 0..100 Scale Ruler", 10, 185, 20, showScaleRuler ? BLACK : LIGHTGRAY);
    DrawText("Segments: ORANGE=Original, RED=Ray, Brown=ROI_Border", 10, 210, 16, GRAY);
    DrawText("Tree edges: DARKBLUE=Original boundary, MAGENTA=Ray boundary", 10, 230, 16, GRAY);
    DrawText(TextFormat("Crossings: %d | Source Face: %d", countCrossingNodes(planarGraph), currentAnalysis.sourceFaceId), 10, 250, 18, MAROON);
    DrawText("Right Click + Drag to Pan | Scroll to Zoom", 10, 275, 20, GRAY);

    Rectangle iterationBox = {(float)width - 250.0f, 10.0f, 180.0f, 32.0f};
    Rectangle roiBoundaryButton = {(float)width - 250.0f, 170.0f, 220.0f, 32.0f};
    const char* roiModeLabel = "neighbors inside";
    switch (roiBoundaryMethod) {
        case ROIBoundaryMethod::NEIGHBORS_INSIDE:
            roiModeLabel = "neighbors inside";
            break;
        case ROIBoundaryMethod::FULL_GRID:
            roiModeLabel = "whole spatial grid";
            break;
        case ROIBoundaryMethod::LOCAL_3X3_AROUND_NODE_CELL:
            roiModeLabel = "local 3x3 around node";
            break;
        default:
            break;
    }

    DrawText("Batch Iterations:", width - 250, 48, 16, DARKGRAY);
    DrawRectangleRec(iterationBox, iterationBoxActive ? Fade(SKYBLUE, 0.2f) : Fade(LIGHTGRAY, 0.35f));
    DrawRectangleLinesEx(iterationBox, 2.0f, iterationBoxActive ? BLUE : DARKGRAY);
    DrawText(iterationInput.empty() ? "0" : iterationInput.c_str(), width - 242, 18, 20, BLACK);
    DrawText("Click box, type digits, press Enter", width - 250, 82, 14, GRAY);
    DrawText(TextFormat("Last batch: %d iterations, %d moves", lastBatchIterations, lastBatchMoves),
             width - 250,
             102,
             14,
             MAROON);
    DrawText(TextFormat("Exploration used: %s", lastBatchExplorationUsed ? "yes" : "no"),
             width - 250,
             120,
             14,
             DARKGRAY);
    if (lastBatchExactCrossings >= 0 && lastBatchNodeBasedCrossings >= 0) {
        DrawText(TextFormat("Batch crossings: exact=%d node=%d diff=%d",
                            lastBatchExactCrossings,
                            lastBatchNodeBasedCrossings,
                            lastBatchCrossingDiff),
                 width - 250,
                 140,
                 14,
                 DARKGREEN);
    }

    DrawText("ROI Boundary Method:", width - 250, 154, 16, DARKGRAY);
    DrawRectangleRec(roiBoundaryButton, Fade(LIGHTGRAY, 0.35f));
    DrawRectangleLinesEx(roiBoundaryButton, 2.0f, DARKGRAY);
    DrawText(TextFormat("Click to switch: %s", roiModeLabel), width - 242, 178, 14, BLACK);
    
    DrawFPS(GetScreenWidth() - 100, 10);
    EndDrawing();
}