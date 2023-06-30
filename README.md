# tree-recovery-by-dynamic-programming
Code for the TPAMI paper: Tree Recovery by Dynamic Programming

## Setup

This project was tested with [Python 3.9.12](https://www.python.org/). Once Python is installed, the required packages can be installed using `pip install -r requirements.txt`. 

## Usage
This code requires two inputs (described in the `Input Data` section):
* Input graph, in the form of a `.ply` file.
* Partition costs, in the form of a `.txt` file.

Use `tree-recovery.py <input graph path> <partition costs path>` to run the optimal Dynamic Programming algorithm. If you wish to use the beam search version, run `tree-recovery.py <input graph path> <partition costs path> <beam width>`

Examples:
* `tree-recovery.py ./datasets/truncated_datasets/rice/graph_01.ply ./datasets/rice/graph_01_costs.txt`
* `tree-recovery.py ./datasets/truncated_datasets/wide/graph_11.ply ./datasets/wide/graph_11_costs.txt 10`

## Input Data

The format is described below. You can also look at the accompanying datasets for examples on the format. 

### Input Graph Format
The input graph is in the form of a [`.ply` file](https://en.wikipedia.org/wiki/PLY_(file_format)). The main requirement is an `edge` element for the edges. This element must have 3 properties: `vertex1` and `vertex2` for the unordered vertices of the edge, and `name`. The `name` is a label for the edge, and will be used to describe the partitions. 

If you plan on using beam search, then you should have a `vertex` element with properties `x` and `y` for coordinates. A `z` property should also be added for 3D datasets.

### Input Costs Format
The input costs file describes both the available partition candidates at each node and their costs.
* The i'th line in the file describes the partitions/costs for the i'th node in the graph.
* Each line has semicolon-separated partition candidates.
* Each partition candidate consists of the following:
  * A space-separated list of integer edge names followed by a comma.
  * For each edge, an integer describing which set it belongs to in the partition. (The order corresponds to the previous list of edge names)
  * For each edge, 0 if it's an outgoing edge, 1 if it's an incoming edge. (The order corresponds to the list of edge names)
  * A floating-point cost.
  * A semicolon, to separated it from the next partition candidate.

For example, suppose the file is as follows:

    2 4, 0 0, 0 0, 0.5;
    4 6 1, 0 0 0, 1 0 0, 0.1; 4 6 1, 0 0 0, 0 1 0, 0.1; 4 6 1, 0 1 1, 1 1 0, 0.3;

Node 0 has one partition candidate. It has the edges 2 and 4. Both edges are assigned to the same set in the partition (the partition is {{2 4}}). Both of them are outgoing edges, meaning this is the root node. The cost of this partition is 0.5. 

Node 1 has 3 partition candidates. They all have the same three edges: 4, 6, and 1.
* For the first candidate, they all belong to the same set (the partition is {{4, 6, 1}}). 4 is an incoming edge, and the other edges outgoing. The cost is 0.1
* For the second candidate, they all belong to the same set. 6 is an incoming edge, and the other edges are outgoing. The cost is also 0.1
* For the third candidate, edge 4 is in one set, and edges 6, 1 are in another set (the partition is {{4}, {6, 1}}). The cost is 0.3

Notes:
* Partitions not listed in this file are *not* considered. If you want a guaranteed solution when running our optimal DP algorithm, make sure to list all possible partitions.
* All cost files must include a "root partition." That is, a node with at least one partition candidate where all edges are outgoing edges.
