#include "logic/Graph.hpp"
#include "logic/PlanarizedGraph.hpp"
#include "logic/SpatialGrid.hpp"
#include "utils/GraphLoader.hpp"
#include "visualization/GraphVisualizer.hpp"
#include "geometry/IntersectionDetector.hpp"
#include <iostream>
#include <vector>

int main() {
    Graph myGraph;

    const char* primaryGraphPath = "src/data/full_dataset/graphswith44nodes/ug5.44";
    const char* fallbackGraphPath = "data/full_dataset/graphswith44nodes/ug5.44";
    bool loaded = GraphLoader::loadFromRomeXML(primaryGraphPath, myGraph) ||
                  GraphLoader::loadFromRomeXML(fallbackGraphPath, myGraph);

    if (loaded) {
        
        auto [minX, minY, maxX, maxY] = myGraph.getBounds(); // Add this helper to Graph class

        // --- STEP 2: TEMPORARY GRID FOR DETECTION ---
        // We use a temporary grid to find the "Egg" (intersections)
        SpatialGrid detectionGrid(minX, maxX, minY, maxY, myGraph.nodes.size());
        
        // Fill grid with original edges
        for (const auto& edge : myGraph.edges) {
            const auto& u = myGraph.nodes[edge.u_id];
            const auto& v = myGraph.nodes[edge.v_id];
            detectionGrid.insertEdge(edge.id, u.x, u.y, v.x, v.y);
        }

        std::vector<PlanarizedGraph::IntersectionData> intersections = IntersectionDetector::findIntersections(myGraph, detectionGrid);

        // --- STEP 3: CONSTRUCT THE FINAL PLANAR GRAPH ---
        // Now the PlanarizedGraph has everything it needs to build the "Chicken"
        PlanarizedGraph pGraph(myGraph, intersections);
        
        // 3. Launch the visualization
        GraphVisualizer visualizer(1280, 720, "Graph Planarization Visualizer");
        visualizer.run(myGraph, pGraph, pGraph.getGrid());
    } else {
        std::cerr << "Failed to load graph." << std::endl;
    }
    
    return 0;
}