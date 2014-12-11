# LSH

[![Build Status](https://travis-ci.org/Keno/LSH.jl.svg?branch=master)](https://travis-ci.org/Keno/LSH.jl)

# Installation

This package requires julia version 0.4 which is currently under development. See http://julialang.org/downloads/ for instructions on how to download a latest nightly release for all major platforms.

After installing julia, this package can be installed with:
```julia
Pkg.clone("https://github.com/Keno/LSH.jl")
# This is currently needed while changes to DataStructures.jl are awaiting upstreaming
Pkg.checkout("DataStructures","kf/lsh")
```

# API

## LSHTable

The main entry point to this package is the `LSHTable` struct, which can be constructed as
```
T = LSHTable(R,hashes,data; progress=true)
```

where `R` is the desired query radius, i.e. the `R` in an `R`-NN data structure, hashes is an array of locality sensitive hash function (see below for how to create these) and `data` is the data set to be added to the locality sensitive hash table. `data` can be specified either as a Vector of vectors, with each contained vector representing a data point, or as a matrix with the coordinates of the points given by the rows of the matrix.

The data structure can be queried via it's `getindex` method:

```
found_points = T[q]
```

where q is the coordinate vector of query points and found_points is the list of found points.

## Hashes

Currently this package supports two kinds of hash functions, both as defined in [AM04]. Both kinds are made up of tuples of hash functions based on `p`-stable distributions, though the first method (dubbed `g`-functions) requires O(kdL) compute whereas the second method (dubbed `u`-functions) can be constructed in O(kdm) where m is appriximately `sqrt(L)`. The hash arrays can be constructed using:

```
hashesg = createHashesAM04(d,w,k,L,R) # Create g functions with the given parameters
hashesu = createUHashesAM04(d,w,k,L,R) # Create u functions with the given parameters
```
