#include "GraphLoader.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <map>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

bool GraphLoader::loadFromRomeXML(const std::string& filepath, Graph& g) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "GraphLoader Error: Could not open file " << filepath << std::endl;
        return false;
    }

    std::string line;
    int currentNode = -1;
    
    // Temporary storage for our two-pass approach
    std::vector<int> parsedNodes;
    std::map<int, std::vector<int>> edgeToNodes; // Map: EdgeID -> List of NodeIDs

    // --- PASS 1: Parse the file ---
    while (std::getline(file, line)) {
        // Reset current node if we hit a closing tag (optional but safe)
        if (line.find("</NODE>") != std::string::npos) {
            currentNode = -1;
            continue;
        }

        // Look for an edge
        if (line.find("<EDGE>") != std::string::npos) {
            std::regex e_regex("<EDGE>\\s*(\\d+)");
            std::smatch match;
            if (std::regex_search(line, match, e_regex)) {
                int edgeId = std::stoi(match[1]);
                if (currentNode != -1) {
                    edgeToNodes[edgeId].push_back(currentNode);
                }
            }
            continue;
        }

        // If we are here, we might be looking at a Node ID line.
        // We look for a line that contains numbers but no XML tags.
        if (currentNode == -1 && line.find('<') == std::string::npos) {
            // Trim leading/trailing whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            if (!line.empty()) {
                try {
                    currentNode = std::stoi(line);
                    parsedNodes.push_back(currentNode);
                } catch (...) {
                    // Not a number, ignore (e.g., empty lines between tags)
                }
            }
        }
    }
    file.close();

    // --- PASS 2: Populate the Graph with a Circle Layout ---
    int numNodes = parsedNodes.size();
    if (numNodes == 0) {
        std::cerr << "GraphLoader Error: No nodes found in file." << std::endl;
        return false;
    }

    // Dynamically scale the circle radius so larger graphs have more room
    double radius = std::max(200.0, numNodes * 15.0); 
    double centerX = radius + 100.0;
    double centerY = radius + 100.0;
    std::map<int, int> nodeIdToIndex;

    for (int i = 0; i < numNodes; ++i) {
        int nId = parsedNodes[i];
        
        // Calculate the angle for this node
        double angle = (2.0 * M_PI * i) / numNodes;
        
        // Calculate coordinates
        double x = centerX + radius * std::cos(angle);
        double y = centerY + radius * std::sin(angle);
        
        g.addNode(nId, x, y);
        nodeIdToIndex[nId] = i;
    }

    // Add the edges
    for (const auto& pair : edgeToNodes) {
        int edgeId = pair.first;
        const auto& connectedNodes = pair.second;
        
        // A valid edge must connect exactly two nodes
        if (connectedNodes.size() == 2) {
            auto uIt = nodeIdToIndex.find(connectedNodes[0]);
            auto vIt = nodeIdToIndex.find(connectedNodes[1]);

            if (uIt != nodeIdToIndex.end() && vIt != nodeIdToIndex.end()) {
                g.addEdge(edgeId, uIt->second, vIt->second);
            } else {
                std::cout << "Warning: Edge " << edgeId << " references unknown node IDs. Skipping.\n";
            }
        } else {
            // Log a warning if the dataset has a malformed edge (e.g., dangling edges)
            std::cout << "Warning: Edge " << edgeId << " does not have exactly 2 endpoints. Skipping.\n";
        }
    }

    return true;
}