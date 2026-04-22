#include "logic/Graph.hpp"
#include "logic/PlanarizedGraph.hpp"
#include "logic/SpatialGrid.hpp"
#include "utils/GraphLoader.hpp"
#include "heuristics/relocation/RelocationHelpers.hpp"
#include "heuristics/RelocationManager.hpp"
#include "geometry/IntersectionDetector.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

namespace {

fs::path resolveSourceRoot() {
    const fs::path current = fs::current_path();
    const std::vector<fs::path> candidates = {
        current,
        current / "..",
        current / ".." / "..",
        current / ".." / ".." / ".."
    };

    std::error_code ec;
    for (const auto& candidate : candidates) {
        const fs::path sourceRoot = candidate / "src";
        if (fs::exists(sourceRoot / "main.cpp", ec) && fs::exists(sourceRoot / "data" / "full_dataset", ec)) {
            return fs::weakly_canonical(sourceRoot, ec);
        }
    }

    return current / "src";
}

std::vector<fs::path> collectDatasetFiles(const fs::path& datasetRoot) {
    std::vector<fs::path> files;
    std::error_code ec;

    auto parseNodeCountFromPath = [](const fs::path& graphPath) {
        for (const auto& part : graphPath) {
            const std::string token = part.string();
            const std::string prefix = "graphswith";
            const std::string suffix = "nodes";

            if (token.rfind(prefix, 0) != 0) continue;
            if (token.size() <= prefix.size() + suffix.size()) continue;
            if (token.substr(token.size() - suffix.size()) != suffix) continue;

            const std::string countText = token.substr(prefix.size(), token.size() - prefix.size() - suffix.size());
            try {
                return std::stoi(countText);
            } catch (...) {
                continue;
            }
        }
        return std::numeric_limits<int>::max();
    };

    auto parseGraphOrdinalFromFilename = [](const fs::path& graphPath) {
        const std::string filename = graphPath.filename().string();
        if (filename.rfind("ug", 0) != 0) return std::numeric_limits<int>::max();

        const size_t dotPos = filename.find('.');
        const std::string ordinalText = filename.substr(2, dotPos == std::string::npos ? std::string::npos : dotPos - 2);
        try {
            return std::stoi(ordinalText);
        } catch (...) {
            return std::numeric_limits<int>::max();
        }
    };

    if (!fs::exists(datasetRoot, ec)) {
        return files;
    }

    for (fs::recursive_directory_iterator it(datasetRoot, ec), end; it != end && !ec; it.increment(ec)) {
        if (!it->is_regular_file(ec)) continue;

        const std::string filename = it->path().filename().string();
        if (filename.rfind("ug", 0) != 0) continue;

        files.push_back(it->path());
    }

    std::sort(files.begin(), files.end(), [&](const fs::path& a, const fs::path& b) {
        const int nodeCountA = parseNodeCountFromPath(a);
        const int nodeCountB = parseNodeCountFromPath(b);
        if (nodeCountA != nodeCountB) return nodeCountA < nodeCountB;

        const int graphOrdinalA = parseGraphOrdinalFromFilename(a);
        const int graphOrdinalB = parseGraphOrdinalFromFilename(b);
        if (graphOrdinalA != graphOrdinalB) return graphOrdinalA < graphOrdinalB;

        return a.generic_string() < b.generic_string();
    });
    return files;
}

bool loadGraphFile(const fs::path& filePath, Graph& graph) {
    graph = Graph{};
    return GraphLoader::loadFromRomeXML(filePath.string(), graph);
}

PlanarizedGraph buildPlanarizedGraph(const Graph& graph) {
    auto [minX, minY, maxX, maxY] = graph.getBounds();
    SpatialGrid detectionGrid(minX, maxX, minY, maxY, graph.nodes.size());

    for (const auto& edge : graph.edges) {
        if (edge.u_id < 0 || edge.v_id < 0 ||
            edge.u_id >= static_cast<int>(graph.nodes.size()) ||
            edge.v_id >= static_cast<int>(graph.nodes.size())) {
            continue;
        }

        const auto& u = graph.nodes[edge.u_id];
        const auto& v = graph.nodes[edge.v_id];
        detectionGrid.insertEdge(edge.id, u.x, u.y, v.x, v.y);
    }

    const auto intersections = IntersectionDetector::findIntersections(graph, detectionGrid);
    return PlanarizedGraph(graph, intersections);
}

int computeExactCrossings(const Graph& graph) {
    if (graph.nodes.empty()) return 0;

    auto [minX, minY, maxX, maxY] = graph.getBounds();
    SpatialGrid detectionGrid(minX, maxX, minY, maxY, graph.nodes.size());

    for (const auto& edge : graph.edges) {
        if (edge.u_id < 0 || edge.v_id < 0 ||
            edge.u_id >= static_cast<int>(graph.nodes.size()) ||
            edge.v_id >= static_cast<int>(graph.nodes.size())) {
            continue;
        }

        const auto& u = graph.nodes[edge.u_id];
        const auto& v = graph.nodes[edge.v_id];
        detectionGrid.insertEdge(edge.id, u.x, u.y, v.x, v.y);
    }

    return static_cast<int>(IntersectionDetector::findIntersections(graph, detectionGrid).size());
}

Graph makeCurrentGraphFromPlanarized(const Graph& originalGraph, const PlanarizedGraph& planarGraph) {
    Graph currentGraph = originalGraph;

    for (auto& node : currentGraph.nodes) {
        if (!planarGraph.hasNode(node.id)) continue;
        const auto& planarNode = planarGraph.getNode(node.id);
        node.x = planarNode.x;
        node.y = planarNode.y;
    }

    return currentGraph;
}

std::string makeCsvRow(const std::string& graphPath,
                       int runIndex,
                       int nodeCount,
                       int edgeCount,
                       int iterations,
                       int initialExactCrossings,
                       int finalExactCrossings,
                       int crossingError,
                       double durationMs) {
    std::ostringstream row;
    row << graphPath << ','
        << runIndex << ','
        << nodeCount << ','
        << edgeCount << ','
        << iterations << ','
        << initialExactCrossings << ','
        << finalExactCrossings << ','
        << crossingError << ','
        << durationMs;
    return row.str();
}

void runBatchDatasetPass() {
    const fs::path sourceRoot = resolveSourceRoot();
    const fs::path datasetRoot = sourceRoot / "data" / "full_dataset";
    const fs::path csvPath = sourceRoot / "batch_results.csv";

    std::ofstream csvFile(csvPath, std::ios::trunc);
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open CSV output file: " << csvPath.string() << '\n';
        return;
    }

    csvFile << "graph_path,run_index,node_count,edge_count,iterations,initial_exact_crossings,final_exact_crossings,crossing_error,duration_ms\n";

    const std::vector<fs::path> graphFiles = collectDatasetFiles(datasetRoot);
    if (graphFiles.empty()) {
        std::cerr << "No dataset graph files found under: " << datasetRoot.string() << '\n';
        return;
    }

    std::error_code ec;
    std::cout << "[Batch] Dataset root: " << datasetRoot.string() << '\n';
    std::cout << "[Batch] CSV output: " << csvPath.string() << '\n';
    std::cout << "[Batch] Graph files found: " << graphFiles.size() << '\n';
    std::cout.flush();

    int graphIndex = 0;
    for (const auto& graphFile : graphFiles) {
        ++graphIndex;
        Graph baseGraph;
        if (!loadGraphFile(graphFile, baseGraph)) {
            std::cout << "[Batch] Skipping unreadable graph: " << graphFile.string() << '\n';
            std::cout.flush();
            continue;
        }

        const int nodeCount = static_cast<int>(baseGraph.nodes.size());
        const int edgeCount = static_cast<int>(baseGraph.edges.size());
        const int iterations = 2 * nodeCount;
        const int initialExactCrossings = computeExactCrossings(baseGraph);
        const fs::path relativeGraphPath = fs::relative(graphFile, datasetRoot, ec);
        const std::string graphPathText = ec ? graphFile.filename().generic_string() : relativeGraphPath.generic_string();

        std::cout << "[Batch] Graph " << graphIndex << "/" << graphFiles.size()
                  << " | file=" << graphPathText
                  << " | n=" << nodeCount
                  << " | m=" << edgeCount
                  << " | initial_exact_crossings=" << initialExactCrossings
                  << " | iterations_per_run=" << iterations
                  << '\n';
        std::cout.flush();

        for (int runIndex = 1; runIndex <= 3; ++runIndex) {
            std::cout << "[Batch] Starting run " << runIndex << "/3 for " << graphPathText << '\n';
            std::cout.flush();

            Graph runGraph = baseGraph;
            PlanarizedGraph planarGraph = buildPlanarizedGraph(runGraph);
            RelocationManager manager(planarGraph);

            const auto start = std::chrono::steady_clock::now();
            for (int step = 0; step < iterations; ++step) {
                manager.performRelocationStep();
            }
            const auto end = std::chrono::steady_clock::now();

            Graph finalGraph = makeCurrentGraphFromPlanarized(runGraph, planarGraph);
            const int finalExactCrossings = computeExactCrossings(finalGraph);
            const int crossingError = finalExactCrossings - initialExactCrossings;
            const double durationMs = std::chrono::duration<double, std::milli>(end - start).count();

            const std::string row = makeCsvRow(graphPathText,
                                               runIndex,
                                               nodeCount,
                                               edgeCount,
                                               iterations,
                                               initialExactCrossings,
                                               finalExactCrossings,
                                               crossingError,
                                               durationMs);

            csvFile << row << '\n';
            csvFile.flush();
            std::cout << row << '\n';
            std::cout.flush();
        }
    }

    std::cout << "[Batch] Completed. Results written to " << csvPath.string() << '\n';
}

} // namespace

int main() {
    runBatchDatasetPass();
    return 0;
}