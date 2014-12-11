using LSH
using Distances
using DataFrames

function constructLSHTable(params,data::Array{Float64,2})
    T = LSHTable(0.6,LSH.createUHashesAM04(params[:d],4.0,params[:k],
                        params[:L],0.6,params[:m]),
    Vector{Float64}[ vec(data[i,:]) for i=1:size(data,1) ]; progress=true)
end

function constructLSHTable(params,data::Vector{Vector{Float64}})
    T = LSHTable(0.6,LSH.createUHashesAM04(params[:d],4.0,params[:k],
                        params[:L],0.6,params[:m]),data; progress=true)
end

function runQueries(T,data,query)
    if typeof(query) == Array{Float64,2}
        query = Vector{Float64}[ vec(query[i,:]) for i=1:size(query,1) ]
    end
    if typeof(data) == Array{Float64,2}
        data = Vector{Float64}[ vec(data[i,:]) for i=1:size(data,1) ]
    end

    # JIT
    r3 = Dict{Any,Any}()
    for (i,q) in enumerate(query)
        local found_points
        time = @elapsed (found_points = T[q])
        nns = length(found_points)
        push!(r3,i,(time,nns))
    end

    r2 = Dict{Any,Any}()
    for (i,q) in enumerate(query)
        local found_points
        time = @elapsed (found_points = T[q])
        nns = length(found_points)
        push!(r2,i,(time,nns))
    end

    exact = Dict{Any,Any}()
    for (i,q) in enumerate(query)
        count = 0
        time = @elapsed for p in data
            if euclidean(p,q) < 0.6
                count += 1
            end
        end
        push!(exact,i,(time,count))
    end

    (r2,exact)
end

function to_dataframes(a,b,c)
    A = DataFrame(pointno=[i for i = 1:length(a)],
            time=[a[i-1][1] for i = 1:length(a)],
            no=[a[i-1][2] for i = 1:length(a)],
            kind=["E2LSH" for i =1:length(a)])
    B = DataFrame(pointno=[i for i = 1:length(b)],
            time=[b[i][1] for i = 1:length(b)],
            no=[b[i][2] for i = 1:length(b)],
            kind=["LSH.jl" for i =1:length(b)])
    C = DataFrame(pointno=[i for i = 1:length(c)],
            time=[c[i][1] for i = 1:length(c)],
            no=[c[i][2] for i = 1:length(c)],
            kind=["Exact" for i =1:length(c)])
    vcat(A,B,C)
end

function compare_to_E2LSH(data, query, jlremote=false)
    if typeof(query) == Array{Float64,2}
        query = Vector{Float64}[ vec(query[i,:]) for i=1:size(query,1) ]
    end
    if typeof(data) == Array{Float64,2}
        data = Vector{Float64}[ vec(data[i,:]) for i=1:size(data,1) ]
    end
    params, e2lshdata = fetch(@spawnat 2 run_e2lsh_on_data(data,query))
    T = constructLSHTable(params,data)
    r2, exact = runQueries(T,data,query)
    e2lshdata, r2, exact
end

function compare_to_E2LSH_mnist(data, query, jlremote=false)
    params, e2lshdata = remotecall_fetch(2,run_mnist_e2lsh_on_data,data,query)
    data = normalized([ trainfeatures(i) for i in data ])
    query = normalized([ testfeatures(i) for i in query ])
    T = constructLSHTable(params,data)
    r2, exact = runQueries(T,data,query)
    e2lshdata, r2, exact
end
