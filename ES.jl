fp_ = pwd()
fn = "RUVrNormalizedCounts.txt"
fp = joinpath(fp_, fn)
# TIMING = 17
# DIM = 128
PB = 0.004
f = CSV.read(fp, DataFrame, delim = '\t', header = 1)
f = f[:, 2:end]

t1 = f[:, 1:3:66]
t2 = f[:, 69:3:110]
t_1 = hcat(t1, t2)
r1 = transpose(Matrix(t_1))

t1 = f[:, 2:3:66]
t2 = f[:, 67:3:110]
t_2 = hcat(t1,t2)
r2 = transpose(Matrix(t_2))


t1 = f[:, 3:3:66]
t2 = f[:, 68:3:110]
t_3 = hcat(t1,t2)
r3 = transpose(Matrix(t_3))

f = CSV.read(fp, DataFrame, delim = '\t', header = 1)
f = f[:, 2:end]

t1 = f[:, 1:3:66]
t2 = f[:, 69:3:110]
t_1 = hcat(t1, t2)
r1 = transpose(Matrix(t_1))
r1 = log.(r1.+1)
ub = maximum(r1)
lb = minimum(r1)
r1N = r1./(ub-lb)

t1 = f[:, 2:3:66]
t2 = f[:, 67:3:110]
t_2 = hcat(t1,t2)
r2 = transpose(Matrix(t_2))
r2 = log.(r2.+1)
ub = maximum(r2)
lb = minimum(r2)
r2N = r2./(ub-lb)

t1 = f[:, 3:3:66]
t2 = f[:, 68:3:110]
t_3 = hcat(t1,t2)
r3 = transpose(Matrix(t_3))
r3 = log.(r3.+1)
ub = maximum(r3)
lb = minimum(r3)
r3N = r3./(ub-lb)

function LargeImpactDetection(AIn, vari)
    A = AIn[1:end-1,:]
    B = AIn[2:end,:]
    X1 = A\B
    var = vari
    flatSol = abs.(vec(X1))
    b = partialsortperm(vec(flatSol), 1:var, rev = true)
    rowN = Array{Int}(undef, var)
    colN = Array{Int}(undef, var)
    for i in 1:var
        tmp = Int(b[i]%4867)
        if tmp == 0
            c = true
        else
            c = false
        end
        rowN[i] = c ? 4867 : tmp
        colN[i] = Int(floor((b[i]-rowN[i])/4867+1))
    end
    rowN, colN
end


mutable struct Ind
    dm:: SparseMatrixCSC{}
    error:: Float64
end
function Base.isless(x::Ind, y::Ind)
    if(x.error < y.error)
        return true
    else
        return false
    end
end



function LRank(num::Int64, p::Float64)
    phi = p
    N = num+1
    a = (2N-phi*(N+1))/(N*(N-1))
    b = (2(phi-1))/((N)*(N-1))
    reso = (N - round((-2a-b+sqrt((2a+b)^2+8*b*rand()))/(2b)))
    while(reso ==0 || reso == N)
        reso = (N - round((-2a-b+sqrt((2a+b)^2+8*b*rand()))/(2b)))
    end
    Int64(reso)
end

function crossover_plain(v1, v2)
    p1 = v1.dm
    p2 = v2.dm
    ccc = p1.&&p2
    c1 = findall(!iszero, p1)
    c2 = findall(!iszero, p2)
    tSN = Int(floor(mean([size(c1,1),size(c2,1)])))
    tt = vcat(c1,c2)
    cc1 = sample(tt, tSN, replace=false)
    cd1 = spzeros(Bool, DIM, DIM)
    cd2 = spzeros(Bool, DIM, DIM)
    cd1[cc1] .=1
    cd2[tt].=1
    cd2[cc1].=0
    dropzeros!(cd2)
    cd1 = cd1.||ccc
    cd2 = cd2.||ccc
    Ind(cd1, 0), Ind(cd2, 0)
end

function crossover_plain__(v1, v2)
    p1 = v1.dm
    p2 = v2.dm
    maxnnz = maximum([nnz(p1), nnz(p2)])
    tmp = (p1 .&& rand(DIM, DIM).<0.4) .|| (p2 .&& rand(DIM, DIM).<0.4)
    c1dm = tmp.|| sprand(Bool, DIM,DIM, PB/10)
    c2dm = tmp.|| sprand(Bool, DIM,DIM, PB/10)
    c1nnz = nnz(c1dm)
    c2nnz = nnz(c2dm)
    if(c1nnz > maxnnz)
        toD = sample(collect(1:c1nnz), c1nnz - maxnnz, replace = false)
        c1dm.nzval[toD] .=0
    end
    if(c2nnz > maxnnz)
        toD = sample(collect(1:c2nnz), c2nnz - maxnnz, replace = false)
        c2dm.nzval[toD] .=0
    end
    dropzeros!(c1dm)
    dropzeros!(c2dm)
    Ind(c1dm, 0), Ind(c2dm, 0)
end

function mutate_flip(pt)
    child = pt.dm
    tm = sprand(Bool, DIM, DIM, 1e-8)
    pt.dm = child.|| tm
end
randspermmat(n) = SparseMatrixCSC(n, n,collect(1:n+1), shuffle(1:n), (1).^rand(Bool,n))
function mutate_permutation(pt)
    if(rand() < 0.01)
        pm = randspermmat(DIM)
        pt.dm = Bool.(pm * pt.dm * pm)
    end
end

function mutate_permutation_ES(pt; mutation_strength = 0.1)
    reso = Ind(Bool.(pt.dm), 0)
    tdn = Int(floor(nnz(reso.dm)*0.1))
    x_ = sample(collect(1:DIM), tdn)
    y_ = sample(collect(1:DIM), tdn)
    tc = findall(!iszero, reso.dm)
    reso.dm[sample(tc, tdn, replace = false)] .= 0
    for (x,y) in zip(x_,y_)
        reso.dm[x,y]=1
    end
    dropzeros!(reso.dm)
    return reso
end



## hyperparameters
# ParentSize = 10
# LRStrength = 0.1
# γ = 25
DATASET_ = r3
r3_mod = DATASET_[1:TIMING, 1:DIM]
# EpochNum_p = 1000
# ParentSize_p = 16
# γ_p = 128
function Run(;out_bool = false, outpath = "", SEED = false, EpochNum = EpochNum_p, runid = -1,ParentSize = ParentSize_p, LRStrength = 0.1, γ = γ_p, DATASET = r3_mod, DIM = DIM)
    ### generate
    flush(stdout)
    println("starting...")
    @printf "μ = %d, λ = %d\n" ParentSize γ*2
    if(SEED == false)
        parents = [Ind(sprand(Bool,DIM, DIM, PB), 0) for _ in 1:ParentSize]
    else
        parents = SEED
    end
    println("evaluating parents...")
    Threads.@threads for i in 1:ParentSize
        evaluate(parents[i], DATASET, off_test)
    end
    ## start iteration
    printBook = Array{Float64}(undef, EpochNum, ParentSize) # record all individuals' fitness per epoch
    nnzbook = Array{Int}(undef, EpochNum) # record number of nonzero numbers per epoch of the best individual
    ## select
    for counter in 1:EpochNum
        @printf "--------Current Epoch is %d------ \n" counter
        sort(parents, rev = true)
        selected_indices=[]
        children = Array{Ind}(undef, γ*2)

        for i in 1:γ
            a = LRank(ParentSize, LRStrength)
            b = LRank(ParentSize, LRStrength)
            while(a==b)
                b = LRank(ParentSize, LRStrength)
            end
            push!(selected_indices, [a,b])
        end
        # println("recombination...")
        Threads.@threads for indx in 1:size(selected_indices,1)
            # println(indx)
            a = selected_indices[indx][1]
            b = selected_indices[indx][2]
            p1 = parents[a]
            p2 = parents[b]
            # children[indx*2-1:indx*2] .= crossover_plain(p1,p2)

            c1, c2 = crossover_plain__(p1,p2)
            c1df = sum(p1.dm .&& c1.dm)+sum(p2.dm .&& c1.dm)
            c2df = sum(p1.dm .&& c2.dm)+sum(p2.dm .&& c2.dm)
            # @printf "c1 sim: %.3f, c2 sim: %.3f \n" c1df c2df
            children[indx*2-1:indx*2] = [c1,c2]
        end
        # println("mutating...")
        mutate_permutation.(children)
        # println("evaluating...")
        Threads.@threads for i in 1:size(children,1)
            # println(i)
            evaluate(children[i], DATASET, off_test)
        end

        total_array = vcat(parents, children)
        sort!(total_array)
        parents = total_array[1:ParentSize]
        hist_book = (map(x->x.error, parents))
        printBook[counter,:] = hist_book
        nnzbook[counter] = nnz(parents[1].dm)
        if(counter%10==0)
            p = plot(1:counter, minimum(printBook[1:counter, :], dims=2),yaxis=:log)
            p = plot!(1:counter, mean(printBook[1:counter, :], dims=2),yaxis=:log)
            p2 = plot(1:counter, nnzbook[1:counter])
            t = plot(p,p2, layout=(2,1), legend = false)
            display(t)
            # @printf "c1 sim: %.3f, c2 sim: %.3f \n" c1df c2df
        end
        @printf "mean is %.4f\n" mean(hist_book)
        @printf "worst is %.4f\n" maximum(hist_book)
        @printf "best is %.4f\n" minimum(hist_book)
        @printf "nnz is %d\n" nnz(parents[1].dm)
    end
    if(out_bool)
        out_filename = joinpath(outpath, string(runid))
        serialize(out_filename*"_hist",[printBook[1:EpochNum, :],nnzbook[1:EpochNum]] ) # record both printbook and nnzbook
        serialize(out_filename*"_models",parents)
    end
    return parents
end


function ES(; SEED = false, EpochNum= 2500, λ = 1024, DATASET = r3_mod, DIM = DIM, PB = 0.004, starting_mutation_strength = 0.1)
    flush(stdout)
    runid = rand(collect(1:99999))
    println("starting...")

    @printf "λ = %d\n" λ
    if(SEED == false)
        parent = Ind(sprand(Bool, DIM, DIM, PB),0)
    else
        parent = SEED
    end
    evaluate(parent, DATASET)
    children = Array{Ind}(undef, λ)
    nnz_hist = Array{Int}(undef, EpochNum)
    performance_hist = Array{Float64}(undef, EpochNum)
    ms = starting_mutation_strength
    unchanged_counter = -1
    for counter in 1:EpochNum
        @printf "--------Current Epoch is %d------ \n" counter
        @printf "ms is %.2f \n" ms
        @printf "unchanged_counter is %d \n" unchanged_counter
        for i in 1:λ
            children[i] = mutate_permutation_ES(parent, mutation_strength = ms)
            evaluate(children[i], DATASET)
        end
        total = vcat(parent, children)
        sort!(total)
        parent = total[1]
        performance_hist[counter] = parent.error
        nnz_hist[counter] = nnz(parent.dm)
        if(counter%10==0)
            p = plot(1:counter, performance_hist[1:counter],yaxis=:log)
            p2 = plot(1:counter, nnz_hist[1:counter])
            t = plot(p,p2, layout=(2,1), legend = false)
            display(t)
            # @printf "c1 sim: %.3f, c2 sim: %.3f \n" c1df c2df
        end
        @printf "error is %.4f\n" performance_hist[counter]
        @printf "nnz is %d\n" nnz_hist[counter]
        if(counter>1)
            if(performance_hist[counter] == performance_hist[counter-1])
                unchanged_counter = unchanged_counter+1
            else
                unchanged_counter = 0
                ms = 0.1
            end
        end
        if(unchanged_counter>ms*ms*1000)
            ms = ms + 0.05
            unchanged_counter = 0
        end
    end
    return parent
end


I_, J_, K_ = findnz(reso)

Z = Array{Bool}(undef, size(K_,1))
Z.=1

inp = sparse(I_ ,J_, Z)
inp_ = spzeros(4867,4867)
d_ = data_test.-off_test
for indx in 1:test_dim
    msk = inp[:,indx].==1
    reso__ = (d_[1:end-1, msk] ) \(d_[2:end, indx] )
    inp_[msk,indx] = reso__
end 
data_test
p = (d_[1:end-1,:])*inp_ .+ off_test
rmsd(data_test[2:end,:], p)


function evaluate(pt, data, offset)
    inp = pt.dm
    inp_ = spzeros(DIM,DIM)
    for indx in 1:DIM
       msk = inp[:,indx].==1
       reso = (data[1:end-1, msk] .- offset[:,msk]) \(data[2:end, indx] .- offset[indx])
       inp_[msk,indx] = reso
    end
    p = (data[1:end-1,:].- off_test)*inp_ .+ off_test
    pt.error = rmsd(data[2:end,:], p)
end


4867*0.0023
0.0023

100 
0.12

test_dim = 4867
vc = 37
timet_test = TIMING
off_test = rand(Uniform(0, 10), 1,test_dim, )
off_test
reso = sprand(test_dim, test_dim,0.0023, n->rand(Uniform(-0.5,0.5),n)) 
reso[diagind(reso)] .= 0
init_test = rand(1, test_dim) 
init_test
data_test = zeros(timet_test, test_dim)
data_test[1,:] = init_test
ubt = off_test.+0.5
lbt = off_test.-0.5

for i in 2:timet_test
    data_test[i, :] = data_test[i-1,:]'*reso 
    # msk = data_test[i,:] .> ubt[:]
    # data_test[i, msk].=ubt[msk]
    # msk = data_test[i,:] .< ubt[:]
    # data_test[i, msk].=lbt[msk]
end


data_test = data_test.+off_test
pl = plot(1:vc, data_test[1:vc, 1:10], legend = false)











pl = plot(1:vc, r3[1:vc, 1:10], legend = false)
pl = plot(1:vc, data_test[1:vc, 1:50], legend = false)


# # data_test[data_test.==0]
# data_test
# # data_test = Array{Float64}(undef, TIMING, DIM)
# # data_test .= 1
p = Run(DATASET =  data_test)


# p[1].dm
# reso

reso
x = findall(!iszero,p[1].dm)
y = findall(!iszero,reso)
ctttttt = 0
for tmp in x
    global  ctttttt
    if(tmp in y)
        ctttttt = ctttttt+1
    end 
end 

println(ctttttt)










































# outdirid_p = string(rand(collect(1:999999)))
# outpath_p = joinpath(pwd(),"data", outdirid_p)
# mkdir(outpath_p)
# println(outpath_p)

# open("data_ids.txt", "a") do file
#     write(file, outdirid_p*"\n")
# end

# ans = "Dim is: " * string(DIM)*"\nTime is: " *string(TIMING) *"\nEpochNum is :" *string(EpochNum_p) *"\nParentSize is :" *string(ParentSize_p) *"\nChildrenSize is :" *string(γ_p*2)

# open(joinpath(outpath_p, "description.txt"), "w") do file
#     write(file, ans)
# end

# for i in 1:10
#     Run(out_bool = true, outpath = outpath_p, EpochNum = EpochNum_p, runid = i )
# end 
