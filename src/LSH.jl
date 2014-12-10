module LSH

####
#
## LSH.jl
#
# This package implements Locality Sensitive Hashing as desribed in [AI06].
# The implementation is guided by the description in [AM04]
#
## Basic Description of the algorithm
#
## Notations and conventions
#
# - Let 𝒫 be the data set we wish to query where each point is of
#   type T𝒫.
#
# - Let H <: HashFunction be the datastructure holding the necessary
#   data to perform a hash of data point yielding a value of type RH.
#   In particular, we need the invariant that
#   ∀ p, typeof(p) = T𝒫 -> H(t)::RH.
#   Implementations for RH === Uint32 are provided in this package.

## References
#
# [AI06] Alexandr Andoni and Piotr Indyk
#   "Near-Optimal Hashing Algorithms for Near Neighbor Problem
#    in High Dimensions"
#
# [AM04] Alexandr Andoni, Mayur Datar, Nicole Immorlica, Piotr Indyk,
#           and Vahab Mirrokni
#   "Locality-Sensitive Hashing Scheme Based on p-Stable Distributions"
#
##

using DataStructures
import DataStructures:  MultiHashDict, AbstractHashDict, HashDict, Unordered, Ordered


abstract HashFunction

export dimension, HashFunction, LSHTable

# ============================= GroupingSet ====================================

import Base: push!, getindex, get
import DataStructures: rehash_item, hashindex

# Groups all elements added to it via `push!(x, k)` into a vector utilizing
# a hash function as the grouping key. All elements inserted into the vector
# so far can be retrieved by looking up any item.

immutable PointGroupingSet{V,GroupingF,RH,HashF,DeriveF,DK,T𝒫} <: AbstractHashDict{V,Vector{V},Unordered}
    group::GroupingF
    𝒫::Vector{T𝒫}
    dict::MultiHashDict{RH,V,Unordered,HashF,DeriveF,DK}
end

function PointGroupingSet{V,GroupingF,RH,HashF,DeriveF,DK,T𝒫}(
        𝒫::Vector{T𝒫},::Type{V},groupF::GroupingF,::Type{RH},hash::HashF,derive::DeriveF,::Type{DK})
    g = PointGroupingSet{V,GroupingF,RH,HashF,DeriveF,DK,T𝒫}(
        groupF,𝒫,MultiHashDict{RH,V,Unordered,HashF,DeriveF,DK}(hash,derive))
    for item = 1:length(𝒫)
        DataStructures._push!(g,g.dict,g.group(g.𝒫[item])=>V(item))
    end
    g
end

rehash(g::PointGroupingSet, sz) = DataStructures._rehash(g,g.dict,sz)

# We can find out what the hash value with this finger print was by passing
# any element of the group back throught the hash
function rehash_item(g::PointGroupingSet, k, item, sz)
    new_key = g.group(g.𝒫[item])
    @assert isequal(k,g.dict.derive(new_key))
    ind = hashindex(g.dict, new_key, sz)
    ind
end

function getindex(g::PointGroupingSet, q)
    key = g.group(q)
    g.dict[key]
end

# ========================= LSH hashindex and derive ===========================

# In the notation of [AM04], we have hashindex == t_1 and derive == t_2

const P = UInt32(UInt64(2)^32 - 5)

immutable ModPHash{RT} <: HashFunction
    # TODO: Make NTuple
    r::Vector{Uint32}
end
dimension(T::ModPHash) = length(T.r)

function check_length(T::ModPHash, N)
    if length(T.r) != N
        error("""Hash function has $(length(T.r)) coefficients, but
                 was called on an input of dimension $N""")
    end
end

function call{RT}(T::ModPHash{RT},z::Vector{Int32})
    check_length(T,length(z))
    result::Uint32 = Uint32(0)
    for i in 1:length(z)
        r = z[i]*UInt32(T.r[i])
        result = ((result + (r % P)) % P) % Uint32
    end
    result % RT
end

function call{RT}(T::ModPHash{RT},z::Int32)
    check_length(T,1)
    (widemul(z,T.r[1]) % P) % RT
end

# ============================== LSHTable ======================================

import Base: push!

# The main data structure on which we may perform query
# contains a collection of k BucketTables
immutable LSHTable{H,RH,DK,V,T𝒫,HashF,DeriveF}
    R::Float64
    𝒫::Vector{T𝒫}
    tables::Vector{PointGroupingSet{V,H,RH,HashF,DeriveF,DK,T𝒫}}
end

# Given a vector of hashes and a dataset, construct an LSHTable
function LSHTable{H,T𝒫,RH,DK}(R::Float64, hashes::Vector{H}, 𝒫::Vector{T𝒫},
        ::Type{RH}, ::Type{DK})

    # Set up the tables.
    tables = [
        PointGroupingSet(
            𝒫,
            Int32,
            h,
            RH,
            ModPHash{Uint32}(rand(Uint32,dimension(h))), # t_1 in [AM04]
            ModPHash{DK}(rand(Uint32,dimension(h))), # t_2 in [AM04]
            DK
        )
        for h in hashes ]

    @show typeof(tables)

    T = LSHTable(R,𝒫,tables)

    T
end

function LSHTable{H,T𝒫}(R::Float64,hashes::Vector{H}, datapoints::Vector{T𝒫})
    LSHTable(R,hashes,datapoints,RH(H),Uint32)
end

using Distances

function getindex{H,RH,DK,V,T𝒫,HashP,DeriveP}(
        T::LSHTable{H,RH,DK,V,T𝒫,HashP,DeriveP},q::T𝒫)
    results = ObjectIdDict()
    tried = falses(length(T.𝒫))
    for table in T.tables
        for p in table[q]
            if !tried[p]
                point = T.𝒫[p]
                if (euclidean(point,q) <= T.R)
                    results[point] = p
                end
                tried[p] = true
            end
        end
    end
    #@show length(tried)
    collect(keys(results))
end

# =============================== LSH Family ===================================
using Distributions

import Distributions: Gaussian, Uniform
import Base: rand

# Implements the LSH family given in [AM04]

immutable AM04HashFamily
    d::Int
    w::Float64
    R::Float64
end

# Represents an h_{a,b} from [AM04]
immutable AM04Hash
    a::Vector{Float64}
    b::Float64
    R::Float64
end

function rand(family::AM04HashFamily)
    #println("Creating new hash")
    AM04Hash(
        rand(Gaussian(0,1.0/family.w),family.d),
        rand(Uniform(0,1)),
        family.R
    )
end

call(h::AM04Hash,v) = floor(Int32,(dot(h.a,v)/h.R + h.b))

RH(::Type{AM04Hash}) = Int32

# The g_i = {h_{1,i}...h_{1,k}}
immutable HashCollection{H}
    hs::Vector{AM04Hash}
end

RH{H}(::Type{HashCollection{H}}) = Vector{H}
dimension(c::HashCollection) = length(c.hs)

call(c::HashCollection,v) = [ h(v) for h in c.hs ]

function createHashesAM04(d,w,k,L,R)
    family = AM04HashFamily(d,w,R)
    [ HashCollection{AM04Hash}([rand(family) for _ = 1:k]) for _ = 1:L ]
end

# Implements the LSH family given in [AI06]

end # module
