using LSH
using Base.Test

# A basic test for LSHTable with a dummy hash function

# To be explicit, we're constructing a test with the following
# parameters to make sure the basic plumbing work (this is not
# supposed to be useful)
#
# - T𝒫 === Int
# - H === DummyHash
# - RH === Uint32
# - DK === Uint32

const P = uint32(uint64(2)^32 - 5)

immutable AdditiveDummyHash <: HashFunction
    z::Int
end
call(T::AdditiveDummyHash,arg) = (T.z+arg)%Int32
LSH.dimension(T::AdditiveDummyHash) = 1

# Store 10 random integers in a table with k=1, L=10
LSHTable([AdditiveDummyHash(z) for z in rand(Int,10)],rand(Int,10),Int32,Int32)
