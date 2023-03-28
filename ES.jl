using CSV
using DataFrames
using Plots
using StatsBase
using StatsFuns
using SparseArrays
using Printf
using Random
using Serialization

#fp = "E:\\julias\\data\\RUVrNormalizedCounts.txt"
fp = "/home/ywang552/research/gene_test/RUVrNormalizedCounts.txt"

TIMING = 17
DIM = 128
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

function evaluate(pt, data)
    inp = pt.dm
    inp_ = spzeros(DIM,DIM)
    for indx in 1:DIM
       msk = inp[:,indx].==1
       reso = data[1:end-1, msk] \ data[2:end, indx]
       inp_[msk,indx] = reso
    end
    p = data[1:end-1,:]*inp_
    pt.error = rmsd(data[2:end,:], p)
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



function Run(;SEED = false, EpochNum = 500, ParentSize = 16, LRStrength = 0.1, γ = 128, DATASET = r3_mod, DIM = DIM)
    ### generate
    flush(stdout)
    runid = rand(collect(1:99999))
    println("starting...")
    @printf "μ = %d, λ = %d\n" ParentSize γ*2
    if(SEED == false)
        parents = [Ind(sprand(Bool,DIM, DIM, PB), 0) for _ in 1:ParentSize]
    else
        parents = SEED
    end
    println("evaluating parents...")
    Threads.@threads for i in 1:ParentSize
        evaluate(parents[i], DATASET)
    end
    ## start iteration
    printBook = Array{Float64}(undef, EpochNum, ParentSize)
    nnzbook = Array{Int}(undef, EpochNum)
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
            evaluate(children[i], DATASET)
        end

        total_array = vcat(parents, children)
        sort!(total_array)
        parents = total_array[1:ParentSize]
        hist_book = (map(x->x.error, parents))
        printBook[counter,:] = hist_book
        nnzbook[counter] = nnz(parents[1].dm)
        # if(counter%10==0)
        #     p = plot(1:counter, minimum(printBook[1:counter, :], dims=2),yaxis=:log)
        #     p = plot!(1:counter, mean(printBook[1:counter, :], dims=2),yaxis=:log)
        #     p2 = plot(1:counter, nnzbook[1:counter])
        #     t = plot(p,p2, layout=(2,1), legend = false)
        #     display(t)
        #     # @printf "c1 sim: %.3f, c2 sim: %.3f \n" c1df c2df
        # end
        @printf "mean is %.4f\n" mean(hist_book)
        @printf "worst is %.4f\n" maximum(hist_book)
        @printf "best is %.4f\n" minimum(hist_book)
        @printf "nnz is %d\n" nnz(parents[1].dm)
        if(counter%100==0)
            # serialize(pwd()*"\\models\\"*string(runid)*"_hist_epoch"*string(counter)*"_("*string(ParentSize)*"+"*string(γ*2)*")",[printBook[1:counter, :],nnzbook[1:counter]] )
            # serialize(pwd()*"\\models\\"*string(runid)*"_models_epoch"*string(counter)*"_("*string(ParentSize)*"+"*string(γ*2)*")",parents)
        end
    end
    return parents
    # serialize(pwd()*"\\models\\"*string(runid)*"_models_epoch"*string(EpochNum)*"_("*string(ParentSize)*"+"*string(γ*2)*")",parents)
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
        # if(counter%10==0)
        #     p = plot(1:counter, performance_hist[1:counter],yaxis=:log)
        #     p2 = plot(1:counter, nnz_hist[1:counter])
        #     t = plot(p,p2, layout=(2,1), legend = false)
        #     display(t)
        #     # @printf "c1 sim: %.3f, c2 sim: %.3f \n" c1df c2df
        # end
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


Run()
# xs = Ind(sprand(Bool, DIM, DIM, PB),0)
# mutate_permutation(xs)
# children = Array{Ind}(undef, 10)
# children[1] = mutate_permutation_ES(xs)


# for i in 1:10
#     Run(EpochNum = 200)
# end
#
# FPFFF = pwd()*"\\models\\"
# P1 = deserialize(FPFFF*"55021_models_epoch200_(16+256)")
# P2 = deserialize(FPFFF*"17611_models_epoch200_(16+256)")
#
#
#
# p1 = P1[1]
# p2 = P2[1]
#
#
#
# r1es = ES(SEED = p1, EpochNum = 200)
# r2es = ES(SEED = p1, EpochNum = 200)
# (r1es.dm.+r2es.dm).==1
# (r1es.dm.+r2es.dm).==2
# end
#
# a, b = crossover_plain__(r1es, r2es)
# evaluate(a, r3_mod)
# evaluate(b, r3_mod)
# P3 = deserialize(FPFFF*"70710_models_epoch200_(16+256)")[1:2]
# P4 = deserialize(FPFFF*"32484_models_epoch200_(16+256)")[1:2]
# P5 = deserialize(FPFFF*"72137_models_epoch200_(16+256)")[1:2]
# P6 = deserialize(FPFFF*"45558_models_epoch200_(16+256)")[1:2]
# P7 = deserialize(FPFFF*"49514_models_epoch200_(16+256)")[1:2]
# P8 = deserialize(FPFFF*"72345_models_epoch200_(16+256)")[1:2]
# P9 = deserialize(FPFFF*"94734_models_epoch200_(16+256)")[1:2]
# P10 = deserialize(FPFFF*"14287_models_epoch200_(16+256)")[1:2]
#
#
# parents_seed = Array{Ind}(undef, 16)
# parents_seed[1:2] = P1
# parents_seed[3:4] = P2
# parents_seed[5:6] = P3
# parents_seed[7:8] = P4
# parents_seed[9:10] = P5
# parents_seed[11:12] = P6
# parents_seed[13:14] = P7
# parents_seed[15:16] = P8
#
# reso = Run(SEED = false, EpochNum = 2000, γ=256, ParentSize = 32)
# reso
# parents_seed[4]== reso[4]
# # #
# P1 = deserialize(FPFFF*"55021_models_epoch200_(16+256)")[1].dm
# P2 = deserialize(FPFFF*"17611_models_epoch200_(16+256)")[1].dm
# P3 = deserialize(FPFFF*"70710_models_epoch200_(16+256)")[1].dm
# P4 = deserialize(FPFFF*"32484_models_epoch200_(16+256)")[1].dm
# P5 = deserialize(FPFFF*"72137_models_epoch200_(16+256)")[1].dm
# P6 = deserialize(FPFFF*"45558_models_epoch200_(16+256)")[1].dm
# P7 = deserialize(FPFFF*"49514_models_epoch200_(16+256)")[1].dm
# P8 = deserialize(FPFFF*"72345_models_epoch200_(16+256)")[1].dm
# P9 = deserialize(FPFFF*"94734_models_epoch200_(16+256)")[1].dm
# P10 = deserialize(FPFFF*"14287_models_epoch200_(16+256)")[1].dm
#
# FPFFF = pwd()*"\\models\\"
# a,b = deserialize(FPFFF*"70710_hist_epoch200_(16+256)")
# plot(1:200, minimum(a,dims=2), yaxis =:log )
# plot!(1:200, mean(a,dims=2), yaxis =:log )
# a = (P1+P2+P3+P4+P5+P6+P7+P8+P9+P10)
# a.==1
# a.==2
# a.==3
