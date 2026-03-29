#pragma once

#include <vector>
#include <string>
#include "logic/Graph.hpp"

class GraphLoader {
public:
    // Parses the GDToolkit XML-like format (.10, .100, etc.)
    // Returns true if successful, false if the file couldn't be read.
    static bool loadFromRomeXML(const std::string& filepath, Graph& g);
};