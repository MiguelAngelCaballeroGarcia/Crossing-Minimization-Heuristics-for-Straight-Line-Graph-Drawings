#include "GraphLoader.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <map>
#include <cmath>
#include <algorithm>

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

    // --- PASS 2: Populate the Graph with a Spiral Layout ---
    int numNodes = parsedNodes.size();
    if (numNodes == 0) {
        std::cerr << "GraphLoader Error: No nodes found in file." << std::endl;
        return false;
    }

    int columns = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(numNodes))));
    int rows = static_cast<int>(std::ceil(static_cast<double>(numNodes) / columns));
    double spacingX = 90.0;
    double spacingY = 90.0;
    double startX = 100.0;
    double startY = 100.0;
    double layoutWidth = columns * spacingX;
    double layoutHeight = rows * spacingY;
    double insetX = spacingX * 0.18;
    double insetY = spacingY * 0.18;
    double centerX = startX + (layoutWidth * 0.5);
    double centerY = startY + (layoutHeight * 0.5);
    double maxRadiusX = std::max(1.0, (layoutWidth * 0.5) - insetX);
    double maxRadiusY = std::max(1.0, (layoutHeight * 0.5) - insetY);
    const double goldenAngle = 2.39996322972865332;
    const double scale = 0.88;
    std::map<int, int> nodeIdToIndex;

    for (int i = 0; i < numNodes; ++i) {
        int nId = parsedNodes[i];

        double normalizedRadius = std::sqrt((static_cast<double>(i) + 0.5) / static_cast<double>(numNodes));
        double angle = i * goldenAngle;

        // Place nodes on a spiral so the layout stays organic and avoids row/column alignment.
        double x = centerX + (std::cos(angle) * maxRadiusX * normalizedRadius * scale);
        double y = centerY + (std::sin(angle) * maxRadiusY * normalizedRadius * scale);
        
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