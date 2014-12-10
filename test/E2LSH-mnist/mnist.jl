# Compare with the sample dataset from E2LSH

using LSH
using Distances
using Base.Test

points = readdlm("mnist1k.dts")
querypoints = readdlm("mnist1k.q")

# d = 784
# w = 4
# Parameters computed by E2LSH was
# k = 14
# L = 153

T = LSHTable(0.6,LSH.createHashesAM04(784,4.0,14,153,0.6),
    Vector{Float64}[ vec(points[i,:]) for i=1:size(points,1) ])

qs = Vector{Float64}[ vec(querypoints[i,:]) for i=1:size(querypoints,1) ]
for (i,q) in enumerate(qs)
    println("=== Querying point $i ===")
    @time found_points = T[q]
    println("Found $(length(found_points)) points. Distances are:")
    for p in found_points
        println("  $(euclidean(p,q))")
    end
end
