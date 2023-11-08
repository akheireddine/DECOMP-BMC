# DECOMP-BMC

This tool aims to combine DeSAT with SAT solvers (namely Painless)
in order to perform Bounded Model Checking (BMC) on DIMACS from BSaLTic encoding (includes LTL spec).

# Dependecies
```
 Read painless README
```
# Install bmc-decomp
make all

# Usage 
  - Sequential context:
    ./bmc-decomp -s=desat -c=1 -decomp=bmc -nleafs=2 file.dimacs
    runs desat with BMC-D decomposition with 2 partitions
    
  - Parallel context:
    ./bmc-decomp -s=minisat-old -c=9 -lbd-limit=4 -bmc-on -decomp=batch -nleafs=5 -shr-strat=1 -d=7 -lbd-limit-itp=4  file.dimacs
    runs a portfolio of 10 threads:
        - 9 minisat1.14 solvers and,
        - 1 random-based decomposition solver using a partition of 5
        - diversification is enable: SparseRandom + native
        - using a simple sharing strategy
        - limit sharing to 4 for both conflict clauses from minisat solvers and LBD<=4 for generated interpolants 