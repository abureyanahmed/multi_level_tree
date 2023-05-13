# A Scalable Method for Readable Tree Layouts

Large tree structures are ubiquitous and real-world relational datasets often have information associated with nodes (e.g., labels or other attributes) and edges (e.g.,weights or distances) that need to be communicated to the viewers. Yet, scalable, easy to read tree layouts are difficult to achieve. We consider tree layouts to be readable if they meet some basic requirements: node labels should not overlap, edges should not cross, edge lengths should be preserved, and the output should be compact. There are many algorithms for drawing trees, although very few take node labels or edge lengths into account, and none optimizes all requirements above. With this in mind, we propose a new scalable method for readable tree layouts. The algorithm guarantees that the layout has no edge crossings and no label overlaps, and optimizes one of the remaining aspects: desired edge lengths and compactness. We evaluate the performance of the new algorithm by comparison with related earlier approaches using several real-world datasets, ranging from a few thousand nodes to hundreds of thousands of nodes. Tree layout algorithms can be used to visualize large general graphs, by extracting a hierarchy of progressively larger trees. We illustrate this functionality by presenting [several map-like visualizations](https://tiga1231.github.io/zmlt/demo/overview.html) generated by the new tree layout algorithm.

## The Readable Tree Layout Algorithm that emphasizes on edge lengths (RTL_L)
See desirable_length_guided/ in this repository

## The Readable Tree Layout Algorithm that emphasizes on companess (RTL_C)
See [[link]](https://github.com/tiga1231/zmlt/) for details

## The Parallel Readable Tree Layout (PRT) Algorithm
See [[link]](https://github.com/khaled-rahman/BatchTree) for details

