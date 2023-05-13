# A Scalable Method for Readable Tree Layouts

Large tree structures are ubiquitous and real-world relational datasets often have information associated with nodes (e.g., labels or other attributes) and edges (e.g.,weights or distances) that need to be communicated to the viewers. Yet, scalable, easy to read tree layouts are difficult to achieve. We consider tree layouts to be readable if they meet some basic requirements: node labels should not overlap, edges should not cross, edge lengths should be preserved, and the output should be compact. There are many algorithms for drawing trees, although very few take node labels or edge lengths into account, and none optimizes all requirements above. With this in mind, we propose a new scalable method for readable tree layouts. The algorithm guarantees that the layout has no edge crossings and no label overlaps, and optimizes one of the remaining aspects: desired edge lengths and compactness. We evaluate the performance of the new algorithm by comparison with related earlier approaches using several real-world datasets, ranging from a few thousand nodes to hundreds of thousands of nodes. Tree layout algorithms can be used to visualize large general graphs, by extracting a hierarchy of progressively larger trees. We illustrate this functionality by presenting [several map-like visualizations](https://tiga1231.github.io/zmlt/demo/overview.html) generated by the new tree layout algorithm.

## The Readable Tree Layout Algorithm that emphasizes on edge lengths (RTL_L)
See desirable_length_guided/ in this repository

To run the code, please open the file named "avoiding_crossing_by_initialization.html" in a browser. The input network needs to be described in a javascript file. A sample input file named "topics_compute_mlst_run_1_postprocessed_5.js" is available in the same directory. The d3 library is available in the folder named "dynamic_update_files". Keep that folder in the working directory and put the javascript file of the d3 library inside that folder. We also need to put some additional files named "area_coverage_mingwei.js", "crossings_initial.js", and "ideal_edge_lenghth_preservation_mingwei.js" in the working directory. These files are also available in the same folder of this repository. These files contain some code to evaluate the algorithms. If we just open the html file, the dataset will be automatically loaded and shown in the browser. Note that if the dataset is large then it can take a few minutes to load. Once loaded, the evaluation of the initial layout will be printed in the console. Then the user can click the start button to run the force-directed algorithm. The user can click the stop button to stop the force-directed algorithm and the evaluation will be printed in the console.

## The Readable Tree Layout Algorithm that emphasizes on companess (RTL_C)
See [[link]](https://github.com/tiga1231/zmlt/) for details

## The Parallel Readable Tree Layout (PRT) Algorithm
See [[link]](https://github.com/khaled-rahman/BatchTree) for details

## Run example code

Download the repository and open the desirable_length_guided/avoiding_crossing_by_initialization.html file in a browser.

## Licence

Copyright 2023 Reyan Ahmed

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
