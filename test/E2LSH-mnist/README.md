## MNIST Test example

This directory contains files and utilities for utilizing the MNIST database of handwritten digits also used in the evaluation of the original E2LSH implementation (and thus serves as the comparison point to this package).

The files in this directory are:

- mnist1k.* - The example dataset included in the E2LSH implementation mnist1k.jl will load this dataset and perform the default queries
- mnist60k.jl - Builds a datastructure using the entire 60k mnist data set, but does not perform any queries (used for experiments)
- compare*.jl - Utilities for reading E2LSH files and running automated comparisons of E2LSH and LSH.jl

The file format of the *.dts and *.q files is dicated by E2LSH and is a space separated list of coordinates, one point per line.
