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
using ProgressMeter
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
    sizehint(g.dict,2*length(𝒫))
    g
end

function push!{V,GroupingF,RH,HashF,DeriveF,DK,T𝒫}(g::PointGroupingSet{V,GroupingF,RH,HashF,DeriveF,DK,T𝒫}, 
        item_no, p=nothing)
    if is(p,nothing)
        DataStructures._push!(g,g.dict,g.group(g.𝒫[item_no])=>V(item_no))
    else
        DataStructures._push!(g,g.dict,g.group(p)=>V(item_no))
    end
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
        r = z[i]*Int64(T.r[i])
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
        ::Type{RH}, ::Type{DK}; progress = false)

    if progress
        p = Progress(length(𝒫), 1)
    end

    @show typeof(𝒫)

    # Set up the tables.
    tables = [
        (pgs = PointGroupingSet(
            𝒫,
            Int32,
            h,
            RH,
            ModPHash{Uint32}(rand(Uint32,dimension(h))), # t_1 in [AM04]
            ModPHash{DK}(rand(Uint32,dimension(h))), # t_2 in [AM04]
            DK
        ); pgs)
        for h in hashes ]

    T = LSHTable(R,𝒫,tables)

    for (i,point) in enumerate(𝒫)
        push!(T,i,point)
        progress && next!(p)
    end

    T
end

function push!(T::LSHTable,i,v)
    v = precompute(T.tables[1].group,v)
    for t in T.tables
        push!(t,i,v)
    end
end

function LSHTable{H,T𝒫}(R::Float64,hashes::Vector{H}, datapoints::Vector{T𝒫}; progress = false)
    LSHTable(R,hashes,datapoints,RH(H),Uint32; progress=progress)
end

using Distances

function getindex{H,RH,DK,V,T𝒫,HashP,DeriveP}(
        T::LSHTable{H,RH,DK,V,T𝒫,HashP,DeriveP},q::T𝒫)
    results = ObjectIdDict()
    tried = falses(length(T.𝒫))
    qp = precompute(T.tables[1].group,q)
    for table in T.tables
        for p in table[qp]
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

# The hashing scheme define in section 3.4.1 of the E2LSH manual

immutable CompositeHashCollection{H}
    first::Uint8
    second::Uint8
    hs::Vector{HashCollection{H}}
end

RH{H}(::Type{CompositeHashCollection{H}}) = Vector{H}
dimension(c::CompositeHashCollection) = 2*dimension(first(c.hs))

immutable PrecomputedHashes{H}
    hs::Vector{HashCollection{H}}
    computed::Vector{Vector{Int32}}
end

call(c::CompositeHashCollection,v) = vcat(c.hs[c.first](v),c.hs[c.second](v))
function call(c::CompositeHashCollection,v::PrecomputedHashes)
    @assert is(c.hs,v.hs)
    vcat(v.computed[c.first],v.computed[c.second])
end

function precompute{H}(hs::Vector{HashCollection{H}},q)
    PrecomputedHashes(hs,[h(q) for h in hs])
end
precompute(h::CompositeHashCollection,q) = precompute(h.hs,q)
precompute(_,q) = q

function createUHashesAM04(d,w,k,L,R,m)
    family = AM04HashFamily(d,w,R)
    hs = [ HashCollection{AM04Hash}([rand(family) for _ = 1:(k/2)]) for _ = 1:m ]
    r = Array(CompositeHashCollection{AM04Hash},div(m*(m-1),2))
    z = 1
    for i = 1:m
        for j = (i+1):m
            r[z] = CompositeHashCollection{AM04Hash}(i,j,hs)
            z += 1
        end
    end
    @assert z-1 == length(r)
    r
end

# Implements the LSH family given in [AI06]

end # module
