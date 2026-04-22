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

## Implementation Status by Phase

### Phase 1: Local Context Extraction (ROI)

**Implemented**
- Neighbor-based ROI bounding from selected node and its neighbors.
- ROI-bounded geometry extraction and clipping.
- ROI border insertion as local boundary segments.
- Ray emission for local subdivision candidates.

**Not Implemented / Partial**
- No separate explicit "cell freezing" subsystem; behavior is effectively achieved through ROI-bounded cell iteration.

### Phase 2: Region Study and Topology Construction

**Implemented**
- Segment splitting at local intersections.
- Face extraction from local arrangement.
- Dual graph construction over inner faces.
- Local analysis pipeline that builds and consumes the dual graph.

**Not Implemented / Partial**
- No major missing item relative to this phase definition.

### Phase 3: Selection and Relocation Update

**Implemented**
- Candidate face selection logic.
- Crossing update plan construction (disappearing/appearing crossing pairs).
- Crossing removal/insertion application.
- Atomic relocation step (move node + reconcile affected crossings).
- Exact per-face global crossing delta evaluation is available.

**Not Implemented / Partial**
- Target selection is still driven primarily by local weight minimization; exact global crossing delta is computed but not yet the primary selection criterion.

### Phase 4: Automated Optimization Loop

**Implemented**
- Batch iteration loop exists in visualization mode and repeatedly executes relocation steps.

**Not Implemented / Partial**
- No simulated annealing acceptance schedule found.
- No explicit convergence policy framework (for example crossing-threshold stop logic or dedicated iteration policy object).
- No standalone non-visual optimizer engine; current iteration loop is UI-driven.

### Phase 5: Observer and Debug Mode

**Implemented**
- ROI overlay rendering.
- Dual-graph debug panel rendering.
- Batch timing and crossing diagnostics shown in UI/console.

**Not Implemented / Partial**
- No structured telemetry pipeline (for example persisted experiment metrics format); current telemetry is mainly UI/console output.

## Build and Run

### Prerequisites
- CMake 3.15+
- C++20-compatible compiler (MSVC, GCC, or Clang)
- `raylib` (resolved via vcpkg or FetchContent)

### Setup & Execute
```powershell
cmake --preset x64-release
cmake --build --preset x64-release
.\build\x64-release\MinimizacionCortes.exe
```