using MNIST

# Extremely ugly quick and dirty
function getE2LSHData(f::IO)
    l1 = "Total time for R-NN query at radius 0.600000 (radius no. 0):\t"
    l2 = "Query point "
    time = 0
    r = Dict{Any,Any}()
    for x in EachLine(f)
        if beginswith(x,l1)
            time = float64(x[endof(l1):end])
        elseif beginswith(x,l2)
            point = int(x[(endof(l2)+1):(search(x,':')-1)])
            if contains(x,"no NNs")
                nn = 0
            else
                nn = int(x[(last(search(x,"found "))+1):(first(search(x," NNs"))-1)])
            end
            push!(r,point,(time,nn))
        end
    end
    r
end

getE2LSHData(filename) = open(getE2LSHData,filename,"r")

function write_e2lsh_points(filename,points::Vector{Vector{Float64}})
    open(filename,"w") do f
        for i = 1:size(points,1)
            for j = 1:size(points[i],1)
                print(f,@sprintf("%0.6lf ",points[i][j]))
            end
            println(f)
        end
    end
end

function write_e2lsh_points(filename,points::Array{Float64,2})
    open(filename,"w") do f
        for i = 1:size(points,1)
            for j = 1:size(points,2)
                print(f,@sprintf("%0.6lf ",points[i,j]))
            end
            println(f)
        end
    end
end

function read_el2sh_parameter_file(filename)
    params = Dict{Symbol,Any}()
    open(filename) do f
        readline(f) # skip the first lne
        while !eof(f)
            x = readline(f)
            if beginswith(x,"R")
                params[:R] = parsefloat(Float64,readline(f))
            elseif beginswith(x,"Success probability")
                params[:P] = parsefloat(Float64,readline(f))
            elseif beginswith(x,"Dimension")
                params[:d] = parseint(readline(f))
            elseif beginswith(x,"Use <u> functions")
                params[:u] = bool(parseint(readline(f)))
            elseif beginswith(x,"k")
                params[:k] = parseint(readline(f))
            elseif beginswith(x,"m")
                params[:m] = parseint(readline(f))
            elseif beginswith(x,"L")
                params[:L] = parseint(readline(f))
            elseif beginswith(x,"R^2") ||
                   beginswith(x,"W") ||
                   beginswith(x,"T") ||
                   beginswith(x,"typeHT")
                readline(f) #skip
            end
        end
        params
    end
end

function run_e2lsh_on_data(data, query; e2lshdir="/home/ubuntu/E2LSH-0.1")
    datafile = tempname()
    write_e2lsh_points(datafile,data)
    queryfile = tempname()
    write_e2lsh_points(queryfile,query)
    paramfile = "$datafile.params"
    local resultdata
    cd(e2lshdir) do
        resultdata = getE2LSHData(`bin/lsh 0.6 $datafile $queryfile`)
    end
    params = read_el2sh_parameter_file(paramfile)
    rm(datafile)
    rm(queryfile)
    rm(paramfile)
    (params,resultdata)
end

function normalized(data)
    [ d / norm(d) for d in data]
end

function run_mnist_e2lsh_on_data(data, query)
    run_e2lsh_on_data(
        normalized([ trainfeatures(i) for i in data ]),
        normalized([ testfeatures(i) for i in query ])
    )
end
