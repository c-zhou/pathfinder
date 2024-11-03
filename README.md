# Pathfinder - genome assembly graph resolution

This is a standalone module extracted from [Oatk](https://github.com/c-zhou/oatk), specifically designed to 
resolve genome assembly graphs by incorporating sequence copy numbers (estimated from sequence coverage) and 
graph structures. It performs an exhaustive search to find a solution that maximizes sequence inclusion while 
adhering to copy number constraints. It is optimized for small assembly graphs, particularly organelle genomes, 
and may produce suboptimal results when applied to larger graphs.

## Installation

You need to have a C compiler, GNU make and zlib development files installed. Download the source code from 
this repo or with `git clone https://github.com/c-zhou/pathfinder.git`. Then type `make` in the source code 
directory to compile.

## Usage

```
Usage: pathfinder [options] <file>[.gfa[.gz]] [<source>[+|-] [<target>[+|-]]]
Options:
  -c INT               maximum copy number of sequences to consider [10]
  -d FLOAT             prefer a circular path if length >= FLOAT * linear length [1.00]
  -p                   do graph partitioning if possible
  -a                   adjust seuqnece copy number estimation by graph structure
  -N INT               maximum number of graph paths to explore [10000000]
 
  --edge-c-tag  STR    edge coverage tag in the GFA file [EC:i] 
  --kmer-c-tag  STR    kmer coverage tag in the GFA file [KC:i] 
  --seq-c-tag   STR    sequence coverage tag in the GFA file [SC:f]
 
  -o FILE              write output to a file [stdout]
  -v INT               verbose level [0]
  --version            show version number

Example: ./pathfinder input.gfa
```

Here are some examples to run.

```
pathfinder TEST/sim_k1001.gfa
pathfinder TEST/sim_k1001.gfa u10
pathfinder TEST/sim_k1001.gfa u10-
pathfinder TEST/sim_k1001.gfa u10- u11+

pathfinder TEST/u2592_f-u2630_f.gfa

pathfinder TEST/u2045_f-u2163_f.gfa
pathfinder -N20000000 TEST/u2045_f-u2163_f.gfa
pathfinder -p TEST/u2045_f-u2163_f.gfa u2045+ u2163+
pathfinder -a -p TEST/u2045_f-u2163_f.gfa u2045+ u2163+

pathfinder -p TEST/u56_r-u144_f.gfa u56- u144
# using two internal nodes u57 and u58
pathfinder -p TEST/u56_r-u144_f.gfa u57- u58-
```

## Major options

`-N` sets the number of solutions to explore in the exhaustive search. While increasing it can sometimes yield 
better results, it is generally not very useful due to the rapid explosion of the solution space.

`-p` may be more useful than `-N` for larger graphs, particularly those with a linear structure and numerous 
entangled local structures. It aims to identify dominators between the source and target nodes, partitioning 
the graph into multiple segments using these dominators as boundaries for individual resolution. To utilize 
this option, you must specify the source and target nodes.

`-a` option adds an additional step to adjust sequence copy numbers by incorporating graph structures. This 
can be beneficial when sequence coverages do not correlate well with copy numbers. However, for very large graphs, 
this option may lead to less accurate estimations.

## Outputs

Below is the output from the command `pathfinder TEST/sim_k1001.gfa`

```
SUBGRAPH node name length coverage copyEstd copyIncl copyDiff
[   1]  u0      10642        198    2    2    0
[   2]  u1      19165         99    1    1    0
[   3]  u2       2714        299    3    3    0
[   4]  u3      22623         99    1    1    0
[   5]  u4      32450         99    1    1    0
[   6]  u5      75279         99    1    1    0
[   7]  u6       7244        199    2    2    0
[   8]  u7      68226         99    1    1    0
[   9]  u8       4187        397    4    4    0
[  10]  u9      47344        100    1    1    0
[  11] u10       2838        199    2    2    0
[  12] u11      51396         99    1    1    0
[  13] u12      15070        100    1    1    0
[  14] u13      51290        100    1    1    0
[  15] u14      65102        100    1    1    0
[  16] u15      11036        200    2    2    0
[  17] u16       1375        200    2    2    0
PATH node name [source=NULL target=NULL nv=27 len=523284 wlen=63725707 circ=ture]
[   1]  u5+
[   2] u16+
[   3]  u9-
[   4]  u8+
[   5] u11+
[   6] u15+
[   7]  u1-
[   8]  u2+
[   9]  u8+
[  10] u10+
[  11] u14+
[  12] u15+
[  13] u12-
[  14] u16-
[  15] u13+
[  16] u10-
[  17]  u8-
[  18]  u2-
[  19]  u0+
[  20]  u4+
[  21]  u6-
[  22]  u3-
[  23]  u0-
[  24]  u2+
[  25]  u8+
[  26]  u7-
[  27]  u6-
```

There are two sections in the output.

The `SUBGRAPH` section shows some basic information for the sequences in the graph, including `name`, `length`, and 
`coverage`, along with the estimated copy numbers (`copyEstd`). It also displays the number of copies of each sequence 
included in the results (`copyIncl`) and the difference from the expected copy numbers (`copyDiff`).

The `PATH` section gives the actual path. The `PATH` header line shows some basic statistics, including `nv` for 
the total number of vertices, `len` for the path length in base pairs, `wlen` for the size weighted by sequence 
coverages, and `circ` to indicate whether the path is circular.
