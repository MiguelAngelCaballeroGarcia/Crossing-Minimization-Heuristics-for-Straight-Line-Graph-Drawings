# MinimizacionCortes

## Project Description

MinimizacionCortes is a C++20 graph-optimization engine focused on minimizing edge crossings in geometric graphs at scale. The implementation utilizes spatial indexing and planarized graph reconstruction to manipulate crossings locally, avoiding expensive global recomputation.

The long-term objective is a fully automated iterative solver for node relocation, using a local topological model derived from the spatial grid.

## Project Status

| Area | Status | Notes |
|---|---|---|
| Data loading | Implemented | Rome-style XML input supported. |
| Intersection detection | Implemented | SpatialGrid-accelerated segment checks. |
| Planarization | Implemented | ORIGINAL/CROSSING node topology maintained. |
| Incremental Spatial Indexing | Implemented | Synchronized DDA-based edge updates. |
| Relocation Logic | **In Progress** | Defining ROI and local geometric extraction. |

## Architecture Overview

- `Graph` (`src/logic`): Base graph representation (nodes and edges).
- `SpatialGrid` (`src/logic`): Geometric acceleration for local edge/node queries.
- `PlanarizedGraph` (`src/logic`): Stitched planarized topology with synchronized grid updates.
- `RelocationManager` (`src/logic`): **(New)** Orchestrates local node movement and ROI calculations.
- `GraphVisualizer` (`src/visualization`): Real-time inspection of graphs and grid layers.

## Algorithm Roadmap

### Phase 1: Local Context Extraction (ROI)
1. **Neighbor-Based Bounding:** Implement logic to identify the relocation **Region of Interest (ROI)** by calculating the minimal square in the `SpatialGrid` covering a node and all its immediate neighbors.
2. **Cell Freezing:** Logic to isolate specific `SpatialGrid` cells within the ROI for localized study.
3. **Geometric Extraction:**
   - Extract all `PlanarEdges` currently occupying the frozen cells ("Static Walls").
   - Implement "Ray Emission": Generate rays from all nodes within the ROI (excluding the variable node) to define potential subdivision boundaries.

### Phase 2: Region Study and Topology Construction
1. **Geometric Structure Analysis:** Analyze intersections between extracted `PlanarEdges` and emitted rays within the ROI boundaries.
2. **Region/Face Construction:** Identify discrete geometric regions formed by the intersection of walls and rays.
3. **Dual Graph Mapping:** - Represent each region as a node in a `DualGraph`.
   - Connect adjacent regions with dual edges, weighted by the crossing cost of the boundary.

### Phase 3: Selection and Relocation Update
1. **Best Region Selection:** Search the `DualGraph` for the region that minimizes the objective function (crossing count).
2. **Atomic Relocation Transaction:**
   - Update physical coordinates of the selected ORIGINAL node.
   - Remove old planar edges/crossings and re-insert new ones via `PlanarizedGraph` helpers.
   - Synchronize `SpatialGrid` indices for the node and its new incident edges.

### Phase 4: Automated Optimization Loop
1. **Iteration Engine:** Wrap the relocation logic in an automated solver that iterates through nodes.
2. **Stochastic Jumping:** Implement a probability-based acceptance logic (Simulated Annealing) to allow moves to "worse" positions, preventing local minima entrapment.
3. **Convergence Policies:** Define exit conditions based on crossing thresholds or iteration caps.

### Phase 5: Observer and Debug Mode (Secondary)
1. **Visual Debugging:** Extend `GraphVisualizer` to render the ROI bounding box, extracted local walls, and candidate regions.
2. **Performance Telemetry:** Track total crossing reduction, iterations per second, and spatial query timings.

## Build and Run

### Prerequisites
- CMake 3.15+
- C++20-compatible compiler (MSVC, GCC, or Clang)
- `raylib` (resolved via vcpkg or FetchContent)

### Setup & Execute
```powershell
cmake --preset x64-debug
cmake --build --preset x64-debug
.\build\x64-debug\MinimizacionCortes.exe