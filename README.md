# MinimizacionCortes

## Project Description

MinimizacionCortes is a C++20 graph-optimization engine focused on minimizing edge crossings in geometric graphs at scale. The current implementation combines spatial indexing and planarized graph reconstruction to detect, represent, and manipulate crossings incrementally, avoiding expensive global recomputation.

The long-term objective is a fully automated iterative solver for node relocation, with deterministic and stochastic search strategies over a locally derived topological model.

## Project Status

| Area | Status | Notes |
|---|---|---|
| Data loading | Implemented | `GraphLoader::loadFromRomeXML(...)` populates `Graph` using Rome-style input. |
| Intersection detection | Implemented | `IntersectionDetector::findIntersections(...)` uses `SpatialGrid` to reduce pair checks. |
| Planarization and stitching | Implemented | `PlanarizedGraph` builds ORIGINAL/CROSSING node topology from intersections. |
| Incremental spatial indexing | Implemented | `SpatialGrid` supports insert/remove for nodes and edges (DDA edge traversal). |
| Interactive visualization baseline | Partial | `GraphVisualizer` supports layers, pan/zoom, and FPS. |
| Optimization engine (automated relocation) | Not implemented | Covered by phased roadmap below. |

## Architecture Overview

- `Graph` (`src/logic`): base graph representation (nodes and edges).
- `GraphLoader` (`src/utils`): dataset parser and graph construction.
- `SpatialGrid` (`src/logic`): geometric acceleration structure for local edge/node queries.
- `IntersectionDetector` (`src/geometry`): robust segment-segment crossing detection.
- `PlanarizedGraph` (`src/logic`): stitched planarized topology with crossing nodes and synchronized grid updates.
- `GraphVisualizer` (`src/visualization`): real-time inspection of original graph, planarized graph, and grid.

## Algorithm Roadmap

The implementation sequence below is the main development path for the solver.

### Phase 1: Iterative Control and Node Selection

1. Implement a dedicated optimization loop (for example, `OptimizationLoop` or `CrossingMinimizer`) for unattended execution.
2. Add pluggable node-selection policies:
   - Round-robin baseline.
   - Random baseline.
   - Highest-crossing-first priority queue.
3. Define convergence policy:
   - No-improvement patience window.
   - Optional hard cap on iterations and/or wall time.

### Phase 2: Candidate Move Space Over SpatialGrid

1. Extract candidate local neighborhoods for a selected ORIGINAL node from `SpatialGrid` occupancy.
2. Subdivide candidate cells into geometric faces induced by existing planar edge segments.
3. Evaluate candidate faces with an explicit objective function:
   - Compute region crossing cost (predicted crossing delta after relocation).
   - Cache local geometry to avoid repeated recomputation.

### Phase 3: Region Connectivity Model

1. Build a `DualGraph` from Phase 2 faces:
   - One dual node per face.
   - One dual edge for each shared face boundary.
2. Assign transition costs to dual edges:
   - Crossing penalty for boundary transitions.
   - Optional extension for geometric regularization terms.

### Phase 4: Stochastic Optimization and Atomic Update

1. Implement target-face search over `DualGraph` (for example shortest path or guided local search).
2. Add stochastic acceptance (Metropolis-Hastings / simulated annealing schedule) to escape local minima.
3. Finalize an atomic relocation transaction for ORIGINAL nodes:
   - Move selected node.
   - Update affected planar sub-edges and crossing nodes.
   - Re-index modified geometry in `SpatialGrid`.
   - Validate topology consistency after each accepted move.

Note: low-level synchronized mutation helpers already exist in `PlanarizedGraph` (`addPlanarEdge`, `removePlanarEdge`, `createCrossing`, `destroyCrossing`), but they are not yet orchestrated by an optimization loop.

### Phase 5: Observer and Debug Mode (Secondary)

1. Extend `GraphVisualizer` with optimization observer features:
   - Current selected node.
   - Candidate faces and chosen target face.
   - Optional `DualGraph` overlay.
2. Add optimization telemetry overlay:
   - Total crossings.
   - Iterations per second.
   - Acceptance ratio.
   - SpatialGrid query/update timings.

## Build and Run

### Prerequisites

- CMake 3.15 or newer.
- C++20-compatible compiler.
- Ninja (used by preset `x64-debug`).
- `raylib` dependency, resolved through vcpkg toolchain or CMake FetchContent fallback.

### Configure

```powershell
cmake --preset x64-debug
```

### Build

```powershell
cmake --build --preset x64-debug
```

### Execute

```powershell
.\build\x64-debug\MinimizacionCortes.exe
```

## Runtime Controls

- `1`: toggle original graph layer.
- `2`: toggle planarized graph layer.
- `G`: toggle spatial grid layer.
- Right mouse drag: pan camera.
- Mouse wheel: zoom camera.

## Dataset Notes

- The entry point currently loads a fixed Rome dataset file from:
  - `src/data/full_dataset/graphswith44nodes/ug5.44` (primary)
  - `data/full_dataset/graphswith44nodes/ug5.44` (fallback)
- For batch experiments, the dataset path should be promoted to a runtime parameter.
