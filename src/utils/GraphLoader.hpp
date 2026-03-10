#pragma once

#include <vector>
#include "Graph.hpp"
#include <string>

// In graph.hpp or a GraphLoader.hpp
struct GraphData {
    Graph graph;
    double minX, maxX, minY, maxY;
};

class GraphLoader {
public:
    static GraphData loadFromDataset(const std::string& filePath) {
        GraphData data;
        // Initialize with extreme values
        data.minX = std::numeric_limits<double>::max();
        data.maxX = std::numeric_limits<double>::lowest();
        data.minY = std::numeric_limits<double>::max();
        data.maxY = std::numeric_limits<double>::lowest();

        // 1. OPEN FILE (pseudo-code depending on your format)
        /* while (reading nodes...) {
            Node n = {id, x, y};
            
            // 2. UPDATE BOUNDS ON THE FLY
            if (x < data.minX) data.minX = x;
            if (x > data.maxX) data.maxX = x;
            if (y < data.minY) data.minY = y;
            if (y > data.maxY) data.maxY = y;

            data.graph.nodes.push_back(n);
        */

        // 3. READ EDGES
        // while (reading edges...) {
        //     data.graph.edges.push_back({id, u, v});
        // }

        return data; // data;
    }
};