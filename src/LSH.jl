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
import DataStructures: AbstractHashDict, HashDict, Unordered, Ordered


abstract HashFunction

export dimension, HashFunction, LSHTable

# ============================= GroupingSet ====================================

import Base: push!

# Groups all elements added to it via `push!(x, k)` into a vector utilizing
# a hash function as the grouping key. All elements inserted into the vector
# so far can be retrieved by looking up any item.

immutable GroupingSet{V,GroupingF,RH,HashF,DeriveF,DK} <: AbstractHashDict{V,Vector{V},Unordered}
    group::GroupingF
    dict::HashDict{RH,Vector{V},Unordered,HashF,DeriveF,DK}
end

function GroupingSet{V,GroupingF,RH,HashF,DeriveF,DK}(::Type{V},groupF::GroupingF,::Type{RH},hash::HashF,derive::DeriveF,::Type{DK})
    GroupingSet{V,GroupingF,RH,HashF,DeriveF,DK}(
        groupF,HashDict{RH,Vector{V},Unordered,HashF,DeriveF,DK}(hash,derive))
end

rehash(g::GroupingSet, sz) = DataStructures._rehash(g,g.dict,sz)

# We can find out what the hash value with this finger print was by passing
# any element of the group back throught the hash
function rehash_item(g::GroupingSet, k, v, sz)
    item = first(v)
    new_key = g.group(item)
    @assert isequeal(g.dict.derive(k),g.dict.derive(new_key))
    hashindex(g.dict, new_key, sz)
end

function push!{T}(g::GroupingSet{T}, v)
    key = g.group(v)

    derived_key = g.dict.derive(key)
    index = DataStructures.ht_keyindex2(g.dict, key, derived_key)

    if index > 0
        arr = g.dict.vals[index]
    else
        arr = Array(T,0)
        DataStructures._setindex!(g.dict, arr, key, derived_key, -index)
    end

    push!(arr,v)
end

# ========================= LSH hashindex and derive ===========================

# In the notation of [AM04], we have hashindex == t_1 and derive == t_2

const P = (uint64(2)^32 - 5)

immutable ModPHash{RT} <: HashFunction
    # TODO: Make NTuple
    r::Vector{Int32}
end
dimension(T::ModPHash) = length(r)

function check_length(T::ModPHash, N)
    if length(T.r) != N
        error("""Hash function has $(length(T.r)) coefficients, but
                 was called on an input of dimension $N""")
    end
end

function call{RT}(T::ModPHash{RT},z::Vector{Int32})
    check_length(T,length(z))
    result::Int = 0
    for i in length(z)
        result = (result + (widemul(z[i],T.r[i]) % P)) % P
    end
    result & RT
end

function call{RT}(T::ModPHash{RT},z::Int32)
    check_length(T,1)
    (widemul(z,T.r[1]) % P) % RT
end

# ============================== LSHTable ======================================

import Base: push!

# The main data structure on which we may perform query
# contains a collection of k BucketTables
immutable LSHTable{H,RH,DK,T𝒫,HashF,DeriveF}
    tables::Vector{GroupingSet{T𝒫,H,RH,HashF,DeriveF,DK}}
end

function LSHTable{H,RH,DK,T𝒫,HashF,DeriveF}(
        tables::Vector{GroupingSet{T𝒫,H,RH,HashF,DeriveF,DK}})
    LSHTable{H,RH,DK,T𝒫,HashF,DeriveF}(tables)
end

# Given a vector of hashes and a dataset, construct an LSHTable
function LSHTable{H,T𝒫,RH,DK}(hashes::Vector{H}, datapoints::Vector{T𝒫},
        ::Type{RH}, ::Type{DK})

    # Set up the tables.
    tables = [
        GroupingSet(
            T𝒫,
            hash,
            RH,
            ModPHash{Int}(rand(Int32,dimension(hash))), # t_1 in [AM04]
            ModPHash{DK}(rand(Int32,dimension(hash))), # t_2 in [AM04]
            DK
        )
        for hash in hashes ]

    T = LSHTable(tables)

    for p in datapoints
        push!(T, p)
    end

    T
end

function push!{H,RH,DK,T𝒫,HashP,DeriveP}(
        T::LSHTable{H,RH,DK,T𝒫,HashP,DeriveP},p::T𝒫)
    for bt in T.tables
        push!(bt,p)
    end
end

end # module
