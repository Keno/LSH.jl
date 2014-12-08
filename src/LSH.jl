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
import DataStructures: AbstractHashDict, Unordered, Ordered


abstract HashFunction

export dimension, HashFunction, LSHTable

# ============================ FixedHashDict ===================================

## Description
# A modified HashDict from DataStructures.jl (which itself is a quadratically
# probed hash table)
# that does not allow rehashing primarily because
# we never need it, but also because we don't actually store enough
# information to do the rehashing.
#
# Furthermore this datastructure allowed the key and the stored key to differ
# related by a derivation function derive(key) => stored_key
#
# In this datastructure function allows varying the hashfunction and derivation
# function in every instance.
#
# The type parameters are:
#
#  - K: typeof(key)
#  - V: typeof(value)
#  - O: Ordering - as usual
#  - DK: typeof(stored_key)
#  - N: Number of slots in the table
#  - DF: The deleter function.
#
## Invariants
#
# derive(K)::DK
#
##

# derivation

immutable DefaultDerivation
end

derive(T::Void) = T
call(d::DefaultDerivation, arg) = derive(arg)

# hashindex function

immutable DefaultHash
end

call(::DefaultHash, arg, n) = hashindex(arg,n)

# Deletion function

immutable CallableIdentity
end

call(c::CallableIdentity, arg) = arg

immutable FixedHashDict{K,V,O<:Union(Ordered,Unordered),DK,N,
        HashF,DeriveF,DeleteF} <: AbstractHashDict{K,V,O}
    slots::Array{Uint8,1}
    keys::Array{DK,1}
    vals::Array{V,1}
    idxs::Array{O,1}
    order::Array{O,1}
    ndel::Int
    count::Int
    hashindex::HashF
    derive::DeriveF
    deleter::DeleteF

    function FixedHashDict(
                hashindex::HashF = DefaultHash(),
                derive::DeriveF = DefaultDerivation(),
                delete::DeleteF = CallableIdentity())
        new(zeros(Uint8,N), Array(DK,N), Array(V,N), Array(O,N),
                Array(O,0), 0, 0, hashindex, derive, delete)
    end
end

# ========================= LSH hashindex and derive ===========================

# In the notation of [AM04], we have hashindex == t_1 and derive == t_2

const P = (uint64(2)^32 - 5)

immutable ModPHash <: HashFunction
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

function call(T::ModPHash,z::Vector{Int32})
    check_length(T,length(z))
    result::Uint32 = 0
    for i in length(z)
        result = ((result + (widemul(z[i],T.r[i]) % P)) % P) % Uint32
    end
end

function call(T::ModPHash,z::Int32)
    check_length(T,1)
    (widemul(z,T.r[1]) % P) % Uint32
end


immutable ModPIndex
    hash::ModPHash
end

call(F::ModPIndex,z,N) = F.hash(z) % N


# ============================= BucketTable ====================================

## Overview
#
# This datastructure associates the table of buckets and the hash function that
# derives the value from the high-dimensional data set
#
## Invariants
#
# H: T𝒫 -> RH
# derive: BucketKey{RH} -> DK
# hashindex(T::BucketKey{RH},N) -> ℤ ⊧ N
##

import Base: push!


immutable BucketTable{H,RH,DK,T𝒫,N,HashP,DeriveP} <:
        Associative{T𝒫, Vector{T𝒫}}
    F::H
    # An unordered hash table {t_2(H(p)) => [p ∈ 𝒫]}
    points::FixedHashDict{RH,Vector{T𝒫},Unordered,DK,N,
                          HashP,DeriveP,CallableIdentity}
end

function BucketTable{H,RH,DK,T𝒫}(F::H,::Type{RH},::Type{DK},::Type{T𝒫},N)
    # Create random hash functions in the family for t_1 and t_2
    idxhash = ModPIndex(ModPHash(rand(Int32,dimension(F))))
    derivehash = ModPHash(rand(Int32,dimension(F)))
    HashT = typeof(idxhash)
    DeriveT = typeof(derivehash)
    BucketTable{H,RH,DK,T𝒫,N,HashT,DeriveT}(F,
        FixedHashDict{RH,Vector{T𝒫},Unordered,DK,N,
            HashT,DeriveT,CallableIdentity}(idxhash,derivehash)
    )
end


function push!{H,RH,DK,T𝒫,N,HashP,DeriveP}(
        BT::BucketTable{H,RH,DK,T𝒫,N,HashP,DeriveP},p::T𝒫)
    # Store H(p) => p in the hash table (with fingerprinting through t_2)
    push!(BT.points,BT.F(p),p)
end

# ============================== LSHTable ======================================

import Base: push!

# The main data structure on which we may perform query
# contains a collection of k BucketTables
immutable LSHTable{H,RH,DK,T𝒫,N,HashP,DeriveP}
    btables::Vector{BucketTable{H,RH,DK,T𝒫,N,HashP,DeriveP}}
end

function LSHTable{H,RH,DK,T𝒫,N,HashP,DeriveP}(btables::Vector{
    BucketTable{H,RH,DK,T𝒫,N,HashP,DeriveP}})
    LSHTable{H,RH,DK,T𝒫,N,HashP,DeriveP}(btables)
end

# Force specialization
immutable dummy{N}; end

function LSHTable{H,T𝒫,RH,DK}(hashes::Vector{H}, datapoints::Vector{T𝒫},
        ::Type{RH}, ::Type{DK})
    N = length(datapoints)
    LSHTable(hashes::Vector{H}, datapoints::Vector{T𝒫},
        RH, DK, dummy{N}())
end

# Given a vector of hashes and a dataset, construct an LSHTable
function LSHTable{H,T𝒫,RH,DK,N}(hashes::Vector{H}, datapoints::Vector{T𝒫},
        ::Type{RH}, ::Type{DK}, ::dummy{N})

    # Set up the BucketTables.
    # As in [AM04], we set the number of buckets
    # in each table to be equal to the number of datapoints in the dataset
    # to guarantee that we can put up with the worst case instance of all
    # datapoints hashing to different indecies.
    btables = [ BucketTable(hash,RH,DK,T𝒫,N)::BucketTable{H,RH,DK,T𝒫,N,
            ModPIndex,ModPHash}
        for hash in hashes ]

    @show typeof(btables)

    T = LSHTable(btables)

    for p in datapoints
        push!(T, p)
    end

    T
end

function push!{H,RH,DK,T𝒫,N,HashP,DeriveP}(
        T::LSHTable{H,RH,DK,T𝒫,N,HashP,DeriveP},p::T𝒫)
    for bt in T.btables
        push!(bt,p)
    end
end

end # module
