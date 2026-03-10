### **README**

**Project Description**

This project implements a high-performance C++ engine designed to minimize edge crossings in large-scale geometric graphs. By utilizing a Spatial Grid and a Local Topological Map, the system identifies and resolves intersections incrementally to ensure maximum scalability. The core objective is to provide a fast, robust framework for topological optimization without the overhead of global recalculations.


**TODO next:**

    1. Implement special case: extremal nodes. Add 4 nodes and 8-12 edges, then don't update the corner nodes, and whenever a extremal node is relocated, update the extremal edges (delete the extremal edges of the updating extremal node, then update its location, then find the new extremal node, then add it the extremal edges).

    2. Run sweep-line algorithm to find crossings coordinates and assign to spatial regions, also find half segments and store them.
