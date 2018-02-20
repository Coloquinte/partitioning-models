# Graph partitioning with a general purpose solver

This is a small script to solve graph partitioning problems with a general purpose solver.

## Principle

Partitioning solvers use a multilevel algorithm: they shrink the size of the graph and apply local search algorithms on the small graph.

This script applies the same principle. However, it does so in a blackbox manner, without looking at the graph's structure. This way, we can tackle complicating constraints and unusual cost functions using a general-purpose solver.

## Running the script

This script was written with LocalSolver in mind, and you will need a LocalSolver license or trial license.
In theory, the script could use any pseudoboolean solver: feel free to adapt it to your solver of choice.


