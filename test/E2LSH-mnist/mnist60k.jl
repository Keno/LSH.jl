using LSH
using Distances
using Base.Test

using MNIST

data = Array(Vector{Float64},60000)
for i = 1:60000
    data[i] = trainfeatures(i)
    data[i] = data[i]/norm(data[i])
end

println("Building data structure")
T = @time LSHTable(0.6,LSH.createUHashesAM04(784,4.0,26,2346,0.6,69),data; progress=true)

nq = 100
query = Array(Vector{Float64},nq)
for i = 1:nq
    query[i] = testfeatures(i)
    query[i] = query[i]/norm(query[i])
end