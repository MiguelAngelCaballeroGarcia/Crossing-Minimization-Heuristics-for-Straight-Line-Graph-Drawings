### **STEP 0 PROJECT**



Define Graph class. ////



Define SpatialGrid class. ////



Iterate over nodes to find bounding box limits and populate Graph object. ////



Special case: extremal nodes. Add 4 nodes and 8-12 edges, then don't update the corner nodes, and whenever a extremal node is relocated, update the extremal edges (delete the extremal edges of the updating extremal node, then update its location, then find the new extremal node, then add it the extremal edges).



Run sweep-line algorithm to find crossings coordinates and assign to spatial regions, also find half segments and store them.



Iterate over nodes to assign spatial regions.



Iterate over half-segments and assign them to regions.

