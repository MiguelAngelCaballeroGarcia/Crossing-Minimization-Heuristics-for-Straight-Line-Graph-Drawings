#pragma once

#include <vector>
#include <cmath>

class Graph;

/**
 * @brief Custom implementation of the Fruchterman-Reingold force-directed algorithm.
 * 
 * This class implements force-directed graph layout using the Fruchterman-Reingold algorithm.
 * It computes repulsive forces between all node pairs and attractive forces along edges.
 * 
 * Supports both:
 * - Exact repulsive force calculation (O(n^2))
 * - Grid-approximated repulsive forces for large graphs (O(n log n))
 */
class FruchtermanReingold {
public:
    /**
     * @struct Vector2D
     * @brief Simple 2D vector for force calculations
     */
    struct Vector2D {
        double x, y;
        
        Vector2D(double x = 0.0, double y = 0.0) : x(x), y(y) {}
        
        Vector2D operator+(const Vector2D& other) const {
            return Vector2D(x + other.x, y + other.y);
        }
        
        Vector2D operator-(const Vector2D& other) const {
            return Vector2D(x - other.x, y - other.y);
        }
        
        Vector2D operator*(double scalar) const {
            return Vector2D(x * scalar, y * scalar);
        }
        
        double dot(const Vector2D& other) const {
            return x * other.x + y * other.y;
        }
        
        double magnitude() const {
            return std::sqrt(x * x + y * y);
        }
        
        Vector2D normalize() const {
            double mag = magnitude();
            if (mag < 1e-10) return Vector2D(0, 0);
            return Vector2D(x / mag, y / mag);
        }
    };

    FruchtermanReingold();
    
    /**
     * @brief Initialize the Fruchterman-Reingold algorithm with drawing area parameters
     * @param boxLength Length of the drawing box
     * @param minX Minimum X coordinate of the drawing area
     * @param minY Minimum Y coordinate of the drawing area
     * @param gridQuotient Grid quotient for approximation (default 2)
     */
    void initialize(double boxLength, double minX, double minY);
    
    /**
     * @brief Calculate exact repulsive forces between all node pairs
     * @param graph The graph structure
     * @param repulsiveForces Output: repulsive force for each node
     */
    void calculateExactRepulsiveForces(const Graph& graph, 
                                      std::vector<Vector2D>& repulsiveForces);
    
    /**
     * @brief Calculate approximated repulsive forces using spatial grid
     * @param graph The graph structure
     * @param spatialGrid The spatial grid for acceleration
     * @param repulsiveForces Output: approximated repulsive force for each node
     */
    void calculateApproximateRepulsiveForces(const Graph& graph,
                                           std::vector<Vector2D>& repulsiveForces);
    
    /**
     * @brief Calculate attractive forces along edges
     * @param graph The graph structure
     * @param attractiveForces Output: attractive force for each node
     */
    void calculateAttractiveForces(const Graph& graph,
                                  std::vector<Vector2D>& attractiveForces);
    
    /**
     * @brief Update node positions based on calculated forces
     * @param graph The graph (positions will be modified)
     * @param repulsiveForces Computed repulsive forces
     * @param attractiveForces Computed attractive forces
     * @param temperature Current temperature (cooling parameter)
     */
    void updateNodePositions(Graph& graph,
                            const std::vector<Vector2D>& repulsiveForces,
                            const std::vector<Vector2D>& attractiveForces,
                            double temperature);
    
    // Parameter setters
    void setOptimalDistance(double distance) { optimalDistance = distance; }
    void setMaxDisplacement(double maxDisp) { maxDisplacement = maxDisp; }
    
    double getOptimalDistance() const { return optimalDistance; }
    double getMaxDisplacement() const { return maxDisplacement; }

private:
    // Drawing area parameters
    double boxLength;
    double minX, minY;
    int gridQuotient;
    int maxGridIndex;
    
    // Algorithm parameters
    double optimalDistance;      // Ideal distance between connected nodes
    double repulsiveStrength;    // Coefficient for repulsive forces
    double attractiveStrength;   // Coefficient for attractive forces
    double maxDisplacement;      // Maximum movement per iteration
    
    /**
     * @brief Calculate repulsive force between two nodes
     * @param pos1 Position of first node
     * @param pos2 Position of second node
     * @return Repulsive force vector (on pos2 away from pos1)
     */
    Vector2D calculateRepulsiveForce(const Vector2D& pos1, const Vector2D& pos2) const;
    
    /**
     * @brief Calculate attractive force between two connected nodes
     * @param pos1 Position of first node
     * @param pos2 Position of second node
     * @return Attractive force vector (on pos2 toward pos1)
     */
    Vector2D calculateAttractiveForce(const Vector2D& pos1, const Vector2D& pos2) const;
    
    /**
     * @brief Calculate forces inside a grid cell (nodes within the same cell)
     * @param cellNodeIndices Indices of nodes in the cell
     * @param graph The graph structure
     * @param repulsiveForces Accumulator for repulsive forces
     */
    void calculateForcesInCell(const std::vector<int>& cellNodeIndices,
                              const Graph& graph,
                              std::vector<Vector2D>& repulsiveForces);
};
