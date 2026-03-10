#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <vector>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;

void printIntersectionCoordinates(const std::vector<Segment>& segments) {
    std::vector<Point> pts;

    // This performs the Bentley-Ottmann sweep-line internally
    CGAL::compute_intersection_points(segments.begin(), segments.end(), std::back_inserter(pts));

    for (const auto& p : pts) {
        // Use CGAL::to_double() to convert exact rational coordinates 
        // to standard doubles for printing/usage.
        std::cout << "Intersection at: (" 
                  << CGAL::to_double(p.x()) << ", " 
                  << CGAL::to_double(p.y()) << ")" << std::endl;
    }
}