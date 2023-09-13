fp_ = pwd()
fn = "RUVrNormalizedCounts.txt"
fp = joinpath(fp_, fn)
f = CSV.read(fp, DataFrame, delim = '\t', header = 1)
f = f[:, 2:end]
p = 0
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

function Ind(x::Ind)
    return Ind(x.dm, x.error)
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
    tm = sprand(Bool, window_size_l, window_size_r, 1e-8)
    pt.dm = xor.(child,  tm)
    return pt 
end
randspermmat(n) = SparseMatrixCSC(n, n,collect(1:n+1), shuffle(1:n), (1).^rand(Bool,n))
function mutate_permutation(pt)
    if(rand() < 0.01)
        pm = randspermmat(DIM)
        pt.dm = Bool.(pm * pt.dm * pm)
    end
end

function mutate_permutation_ES(pt, mutation_strength = 0.1)
    reso = Ind(Bool.(pt.dm), 0)
    tdn = Int(floor(nnz(reso.dm)*mutation_strength))
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


# function fastmut(pt, ms, ub, lb)
#     local nr
#     cond = true
#     while(cond)
#         a = pt.dm .&& rand(DIM, DIM).<(1-ms)
#         nr = dropzeros(a.|| rand(DIM,DIM).<PB*ms/((1-PB)+PB*ms))
#         if(nnz(nr) < ub && nnz(nr)> lb)
#             cond = false
#         end 
#     end 
#     return Ind(nr,0)
# end 

function fastmut(pt, ms)
    local nr
    a = pt.dm .&& rand(DIM, DIM).<(1-ms)
    nr = dropzeros(a.|| rand(DIM,DIM).<PB*ms/((1-PB)+PB*ms))
    return Ind(nr,0)
end 

function fastmut_window(pt, ms)
    it, jt, kt, = findnz(pt.dm)
    tt = nnz(pt.dm)
    ts = Int(floor(ms*tt))
    sp = sample(1:tt, tt - ts, replace=false)
    for toreplace in sp 
        l = rand(1:window_size_l)
        r = rand(1:window_size_r)
        while((l,r) in zip(it, jt))

            l = rand(1:window_size_l)
            r = rand(1:window_size_r)

        end 
        it[toreplace] = l 
        jt[toreplace] = r
    end 
    rs = sparse(it,jt, kt, window_size_l, window_size_r)
    # rs[diagind(rs)] .= 0
    # dropzeros!(rs)
    return Ind(rs, 0)
end 


function imba_mutation(ptr)
    mut = 0.1
    a = zeros(size(ptr.dm.nzval))
    b = similar(ptr.dm.nzval, Vector{Int16})
    idx = 1
    for (x,y,z) in zip(findnz(ptr.dm)...)
        local pred
        r = data_test[1:end-1, x] \ data_test[2:end, y]
        pred = data_test[1:end-1, x] * r .+ off_test[y]
        a[idx] = rmsd(pred, data_test[2:end, y])
        b[idx] = [x,y]
        idx = idx + 1
    end 
    perm = sortperm(a)
    b .= b[perm]
    c = zeros(Int, size(b,1), 2)
    c[1:end, 1] .= map(x->x[1], b)
    c[1:end, 2] .= map(x->x[2], b)
    tar = Int(floor(size(b,1)*mut))

    nr = spzeros(Bool, DIM,DIM)
    l = zeros(Int, size(b))
    r = zeros(Int, size(b))
    l[1:tar] = c[1:tar,1]
    r[1:tar] = c[1:tar,2]
    for i in eachindex(b)
        if(i>tar)
            a_,b_ = rand(1:DIM), rand(1:DIM)
            while((a_,b_) in zip(l,r))
                a_,b_ = rand(1:DIM), rand(1:DIM)
            end
            l[i] = a_
            r[i] = b_
        end 
    end 

    for (x,y) in zip(l,r)
        nr[x,y] = 1
    end 
    Ind(nr, 0)
end 

function mut_mut(ptr)
    p = ptr.dm
    nr = spzeros(Bool, DIM, DIM)
    mut = 0.1
    colCounter = 1
    for col in eachcol(p)
        tp = Int(sum(col))
        s = Int(ceil(tp*mut))
        ta = Array{Float64}(undef, tp)
        ra = Array{Int16}(undef, tp)
        ca = ones(Int16,tp)*colCounter
        rowCounter = 1
        for index in eachindex(col)
            if(p[index,colCounter] == 1)
                r = data_test[1:end-1, index] \ data_test[2:end, colCounter]
                pred = data_test[1:end-1, index] * r .+ off_test[colCounter]
                ta[rowCounter] = rmsd(pred, data_test[2:end, colCounter])
                ra[rowCounter] = index
                rowCounter = rowCounter + 1
            end 
        end 
        perm = sortperm(ta)
        ra = ra[perm]
        for index in s+1:tp
            newRow = rand(1:DIM)
            while(newRow in ra)
                newRow = rand(1:DIM)
            end 
            ra[index] = newRow
        end 
        for (x,y) in zip(ra,ca)
            nr[x,y] = 1
        end 
        colCounter = colCounter + 1
    end 

    return Ind(nr,0) 
end 

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
            p3 = plotSol(parents[1])

            t = plot(p,p2, p3, layout=(3,1), legend = false)
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


function ES(; SEED = false, EpochNum= 2500, λ = 1024, DATASET = r3_mod, DIM = DIM, PB = PB, starting_mutation_strength = 0.8)
    flush(stdout)
    runid = rand(collect(1:99999))
    println("starting...")

    @printf "λ = %d\n" λ
    if(SEED == false)
        parent = Ind(sprand(Bool, DIM, DIM, PB),0)
    else
        parent = SEED
    end
    evaluate(parent, DATASET, off_test)
    children = Array{Ind}(undef, λ)
    nnz_hist = Array{Int}(undef, EpochNum)
    performance_hist = Array{Float64}(undef, EpochNum)
    ms = starting_mutation_strength
    ub = Int(floor(nnz(parent.dm)*1.008))
    lb = Int(floor(nnz(parent.dm)*0.992))
    unchanged_counter = -1
    for counter in 1:EpochNum
        @printf "--------Current Epoch is %d------ \n" counter
        @printf "ms is %.2f \n" ms
        @printf "unchanged_counter is %d \n" unchanged_counter
        for i in 1:λ
            if(rand()<0.5)
                children[i] = fastmut(parent, ms, ub, lb)
            else
                children[i] = mut_mut(parent)
                children[i] = imba_mutation(children[i])
            end 
            evaluate(children[i], DATASET,off_test)
        end
        total = vcat(parent, children)
        # total = vcat(children)
        sort!(total)

        parent = total[1]
        performance_hist[counter] = parent.error
        nnz_hist[counter] = nnz(parent.dm)
        # if(counter < 20)
        #     t = plot(plotSol(parent), title = "epoch number "*string(counter))
        #     display(t)

        # elseif(counter%10==0)
        #     # p = plot(1:counter, performance_hist[1:counter],yaxis=:log)
        #     # p2 = plot(1:counter, nnz_hist[1:counter])
        #     # hline!([ub ub])
        #     # hline!([lb lb])
        #     # t = plot(p,p2, layout=(2,1), legend = false)
        #     # display(t)
        #     t = plot(plotSol(parent), title = "epoch number "*string(counter))
        #     display(t)
        #     # @printf "c1 sim: %.3f, c2 sim: %.3f \n" c1df c2df
        # end
        if(unchanged_counter == 0 || unchanged_counter == -1)
            t = plot(plotSol(parent), title = "epoch number "*string(counter))
            display(t)
        end 


        @printf "error is %.4f\n" performance_hist[counter]
        @printf "nnz is %d\n" nnz_hist[counter]
        if(counter>1)
            if(performance_hist[counter] == performance_hist[counter-1])
                unchanged_counter = unchanged_counter+1
            else
                unchanged_counter = 0
                ms = starting_mutation_strength
            end
        end
        if(unchanged_counter>ms*ms*1000)
            ms = ms + 0.05
            unchanged_counter = 0
        end
    end
    return parent
end


function evaluate(pt, data, offset)
    inp = pt.dm
    inp_ = spzeros(DIM,DIM)
    data_ = data.-offset
    for indx in 1:DIM
       msk = inp[:,indx].==1
       msk = CartesianIndices(msk)[msk]
       reso_ = (data_[1:end-1, msk] ) \ (data_[2:end, indx] )
       inp_[msk,indx] = reso_
    end
    p = (data_[1:end-1,:])*inp_ .+ offset
    pt.error = rmsd(data[2:end,:], p)
end


function calculate(pt, data, offset)
    inp = pt.dm
    inp_ = spzeros(DIM,DIM)
    data = data.-offset
    for indx in 1:DIM
       msk = inp[:,indx].==1
       msk = CartesianIndices(msk)[msk]
       reso_ = (data[1:end-1, msk] ) \ (data[2:end, indx] )
       inp_[msk,indx] = reso_
    end
    p = (data[1:end-1,:])*inp_ .+ offset
    return p
end

function evaluate_gene(pt, data, offset, i)
    inp = pt.dm
    inp_ = spzeros(DIM,DIM)
    data_ = data.-offset
    for indx in 1:DIM
       msk = inp[:,indx].==1
       msk = CartesianIndices(msk)[msk]
       reso_ = (data_[1:end-1, msk] ) \ (data_[2:end, indx] )
       inp_[msk,indx] = reso_
    end
    p = (data_[1:end-1,:])*inp_ .+ offset
    return rmsd(data[2:end,i], p[:,i])
end
function plotSol(ptr)
    m = ptr.dm
    c1 = 0
    c2 = 0 
    c3 = 0
    a = map(x->RGB(0,0,0), m)
    for i in eachindex(r_)
        if(m[i] == 1 && r_[i] == 1)
            a[i]=RGB(0,1,0)
            c1 = c1+1
        end 
        if(m[i] == 1 && r_[i] == 0)
            a[i]=RGB(1,0,0)
            c2 = c2+1
        end 
        if(m[i] == 0 && r_[i] == 1 )
            a[i] = RGB(0.571945,0.42736,0.548254)
            c3 = c3+1
        end 
    end 
    t1 = "c1:" *string(c1)
    t2 = "c2:" *string(c2)
    t3 = "c3:" *string(c3)
    res = plot(a, axis = nothing, annotations = (DIM, 0, Plots.text(t1, :left)))
    annotate!(DIM, 40, Plots.text(t2, :left))
    annotate!(DIM, 80, Plots.text(t3, :left))
    return res
end 

function plotSol_window(ptr, is, js, windowsize_l, windowsize_r)
    m = ptr.dm

    c1 = 0
    c2 = 0 
    c3 = 0
    r = r_[1+is*windowsize_l:(is+1)*windowsize_l, 1+js*windowsize_r: (js+1)*windowsize_r]
    a = map(x->RGB(0,0,0), m)
    for i in eachindex(r)
        if(m[i] == 1 && r[i] == 1)
            # println("zz", evaluate_singleGene(i, data_test, off_test))

            a[i]=RGB(0,1,0)
            c1 = c1+1
        end 
        if(m[i] == 1 && r[i] == 0)
            if(i[1] != i[2])
                a[i]=RGB(1,0,0)
            end
            # println(evaluate_singleGene(i, data_test, off_test))
            c2 = c2+1
        end 
        if(m[i] == 0 && r[i] == 1 )
            a[i] = RGB(0.571945,0.42736,0.548254)
            c3 = c3+1
        end 
    end 
    t1 = "c1:" *string(c1)
    t2 = "c2:" *string(c2)
    t3 = "c3:" *string(c3)
    res = plot(a, axis = nothing, annotations = (windowsize_r, 1, Plots.text(t1, :left)))
    annotate!(windowsize_r, 0.3 * window_size_l, Plots.text(t2, :left))
    annotate!(windowsize_r, 0.6 * window_size_l, Plots.text(t3, :left))
    return res
end 


function calSol_window(ptr, is, js, windowsize_l, windowsize_r)
    m = ptr.dm

    c1 = 0
    c2 = 0 
    c3 = 0
    r = r_[1+is*windowsize_l:(is+1)*windowsize_l, 1+js*windowsize_r: (js+1)*windowsize_r]
    a = map(x->RGB(0,0,0), m)
    for i in eachindex(r)
        if(m[i] == 1 && r[i] == 1)
   
            a[i]=RGB(0,1,0)
            c1 = c1+1

        end 
        if(m[i] == 1 && r[i] == 0)
            if(i[1] != i[2])
                a[i]=RGB(1,0,0)
                c2 = c2+1
            end
        end 
        if(m[i] == 0 && r[i] == 1 )
            a[i] = RGB(0.571945,0.42736,0.548254)
            c3 = c3+1
        end 
    end 
    return c1, c2
end 


function green_window(ptr, is, js, windowsize_l, windowsize_r)
    m = ptr.dm

    c1 = 0
    r = r_[1+is*windowsize_l:(is+1)*windowsize_l, 1+js*windowsize_r: (js+1)*windowsize_r]
    for i in eachindex(r)
        if(m[i] == 1 && r[i] == 1)
            c1 = c1+1
        end 
    end 
    return c1
end 


function evaluate_sw(pt, data, offset, is, js, window_size_l, window_size_r)
    inp = pt.dm
    inp_ = spzeros(window_size_l,window_size_r)
    data_ = data.-offset
    for indx in 1:window_size_r
        msk = inp[:,indx] .== 1
        msk = CartesianIndices(msk)[msk]
        msk_ = similar(msk)
        for i in eachindex(msk)
            msk_[i] = msk[i] + CartesianIndex(window_size_l*is)
        end 
        reso = data_[1:end-1, msk_] \ data_[2:end, indx+window_size_r*js] 
        inp_[msk,indx] = reso
    end 
    p = (data_[1:end-1, 1+window_size_l*is : window_size_l*(is+1)])*inp_ .+offset[1+window_size_r*js:window_size_r*(js+1)]'
    pt.error = rmsd(data[2:end,1+window_size_r*js:window_size_r*(js+1)], p)
end 

function evaluate_singleGene(indx, data, offset)
    data_ = data.-offset
    rid, cid = indx[1], indx[2]
    rs = data_[1:end-1,rid] \ data_[2:end, cid]
    rmsd(data_[1:end-1,rid]*rs.+offset[cid]', data[2:end, cid])
end 


function evaluate_sw_printf(pt, data, offset, is, js, window_size_l, window_size_r)
    inp = pt.dm
    inp_ = spzeros(window_size_l,window_size_r)
    data_ = data.-offset
    for indx in 1:window_size_r
        msk = inp[:,indx] .== 1
        msk = CartesianIndices(msk)[msk]
        msk_ = similar(msk)
        for i in eachindex(msk)
            msk_[i] = msk[i] + CartesianIndex(window_size_l*is)
        end 
        reso = data_[1:end-1, msk_] \ data_[2:end, indx+window_size_r*js] 
        inp_[msk,indx] = reso
    end 
    p = (data_[1:end-1, 1+window_size_l*is : window_size_l*(is+1)])*inp_ .+offset[1+window_size_r*js:window_size_r*(js+1)]'
    return rmsd(data[2:end,1+window_size_r*js:window_size_r*(js+1)], p)
end 

function evaluate_cheat(ptr, is, js, windowsize_l, windowsize_r)
    r = r_[1+is*windowsize_l:(is+1)*windowsize_l, 1+js*windowsize_r: (js+1)*windowsize_r]
    
    ptr.error = sum(xor.(ptr.dm, r) )

end 

function evaluate_cheat_print(ptr, is, js, windowsize_l, windowsize_r)
    r = r_[1+is*windowsize_l:(is+1)*windowsize_l, 1+js*windowsize_r: (js+1)*windowsize_r]
    
    return sum(xor.(ptr.dm, r) )

end 

function ES_sliding(; SEED = false, is = 0, js = 0, EpochNum= 2500, λ = 1024, DATASET = r3_mod, DIM_ = DIM, PB_ = PB, starting_mutation_strength = 0.8, windsize_l = 32, windsize_r = DIM)
    flush(stdout)
    println("starting...")
    @printf "λ = %d\n" λ
    if(SEED == false)
        parent = Ind(sprand(Bool, windsize_l, windsize_r, PB),0)
    else
        parent = SEED
    end
    # evaluate_sw(parent, DATASET, off_test, is, js, windsize_l, windsize_r)
    evaluate_cheat(parent, is, js, windsize_l, windsize_r)    

    children = Array{Ind}(undef, λ)
    nnz_hist = Array{Int}(undef, EpochNum)
    performance_hist = Array{Float64}(undef, EpochNum)
    ms = starting_mutation_strength
    ub = Int(floor(nnz(parent.dm)*1.008))
    lb = Int(floor(nnz(parent.dm)*0.992))
    unchanged_counter = -1
    for counter in 1:EpochNum
        @printf "--------Current Epoch is %d------ \n" counter
        @printf "ms is %.2f \n" ms
        @printf "unchanged_counter is %d \n" unchanged_counter 
        for i in 1:λ
            # children[i] = fastmut_window(parent, ms)
            children[i] = mutate_flip(parent)
            # evaluate_sw(children[i], DATASET, off_test,is, js, windsize_l, windsize_r)
            evaluate_cheat(children[i], is, js, windsize_l, windsize_r)    
        end
        total = vcat(parent, children)
        # total = vcat(children)
        sort!(total)

        parent = total[1]
        performance_hist[counter] = parent.error
        nnz_hist[counter] = nnz(parent.dm)

        if(unchanged_counter == 0 || unchanged_counter == -1)
            t = plot(plotSol_window(parent, is, js, windsize_l, windsize_r), title = "epoch number "*string(counter))
            display(t)
        end 


        @printf "error is %.8f\n" performance_hist[counter]
        @printf "nnz is %d\n" nnz_hist[counter]
        if(counter>1)
            if(performance_hist[counter] == performance_hist[counter-1])
                unchanged_counter = unchanged_counter+1
            else
                unchanged_counter = 0
                ms = starting_mutation_strength
            end
        end
        if(unchanged_counter>ms*ms*1000)
            ms = ms + 0.05
            unchanged_counter = 0
        end
    end
    return parent
end 

function distance(x, y, zs)
    t = sum(x.dm .|| y.dm) - sum(x.dm .&& y.dm)
    if(t < zs)
        return 1 - ((t)/zs)^2
    else
        return 0 
    end 
end

p1 = Ind(sprand(Bool, window_size_l, window_size_r, PB),0)
p2 = Ind(sprand(Bool, window_size_l, window_size_r, PB),0)

function crossover_withMutation(ptr1, ptr2, windl, windr, inheretRate)
    r = spzeros(Bool,windl, windr)
    genePool_p1 = findall(!iszero, ptr1.dm)
    genePool_p2 = findall(!iszero, ptr2.dm)
    genePool_total = vcat(genePool_p1, genePool_p2)
    geneNum_pL = minimum([nnz(ptr1.dm), nnz(ptr2.dm)])
    geneNum_pR = maximum([nnz(ptr1.dm), nnz(ptr2.dm)])
    geneNum_childTotal = sample(geneNum_pL:geneNum_pR)
    geneNum_childInhereted = Int(floor(geneNum_childTotal*inheretRate))
    genePool_childInhereted = unique(sample(genePool_total, geneNum_childInhereted, replace = false))
    tmp = size(genePool_childInhereted,1)
    if(tmp < geneNum_childInhereted)
        for i in 1:10
            genePool_childInhereted = unique(vcat(genePool_childInhereted, sample(genePool_total, geneNum_childInhereted - tmp)))
            tmp = size(genePool_childInhereted, 1)
            if(tmp >= geneNum_childInhereted)
                break
            end 
        end 
    end 
    r[genePool_childInhereted] .= 1

    for i in 1:geneNum_childTotal-geneNum_childInhereted
        x = sample(1:windl)
        y = sample(1:windr)
        c = 1
        while((r[x,y] == 1 ||  x==y) || c < 100)
            x = sample(1:windl)
            y = sample(1:windr)
            c = c+1
        end 
        r[x,y] = 1
    end
    dropzeros!(r)
    return Ind(r, 0)
end 

function crossover_deterministic(ptr1, ptr2, windl, windr, inheretRate)
    r = ptr1.dm .&& ptr2.dm
    r2 = xor.(ptr1.dm, ptr2.dm)

    mini = minimum([nnz(ptr1.dm), nnz(ptr2.dm)])
    maxi = maximum([nnz(ptr1.dm), nnz(ptr2.dm)])
    geneNum_children = sample(mini:maxi) - nnz(r)
    gene_setTotal = findall(!iszero, r2)
    shuffle!(gene_setTotal)
    msk = gene_setTotal[1:geneNum_children]
    r[msk].=1

    dropzeros!(r)
    return Ind(r,0)
end 


function crossover_deterministic_(ptr1, ptr2, windl, windr, inheretRate)
    r = ptr1.dm .&& ptr2.dm
    r2 = xor.(ptr1.dm, ptr2.dm)
    if(sum(r2) !=0)
        common_num = nnz(r)
        toIgnore = Int(floor(common_num*0.05))
        common_set = findall(!iszero, r)
        msk = sample(common_set, toIgnore, replace = false)
        r[msk].=0
        dropzeros!(r)
    end 

    mini = minimum([nnz(ptr1.dm), nnz(ptr2.dm)])
    maxi = maximum([nnz(ptr1.dm), nnz(ptr2.dm)])
    gene_setTotal = findall(!iszero, r2)
    geneNum_children = minimum([maximum([sample(mini:maxi) - nnz(r), 0]), size(gene_setTotal,1)])
    shuffle!(gene_setTotal)
    msk = gene_setTotal[1:geneNum_children]
    r[msk].=1

    dropzeros!(r)
    return Ind(r,0)
end 


function ES_cheat(; verbose = verbose_, seeded = false, seed = nothing,  is = is_, js = js_, epochNum = EpochNum_p, μ = μ_, λ = λ_, DATASET = data_test, PB = PB_, windrow = window_size_l, windcol = window_size_r, p = p_, ms = ms_, inheretRate = inheretRate_, eli_num = eli_num_, restart_num = restart_num_)
    zone_size = windrow*windcol*PB*0.2
    if(verbose)
        @printf "starting...\n"
        @printf "(μ+λ) is (%d + %d) \n" μ  λ 
        @printf "the ms is %.2f\n" ms
        @printf "seed is %s\n" seeded
    end 
    rv = true     
    if(seeded)
        parents = seed
    else
        parents = Array{Ind}(undef, μ)
        for i in 1:μ
            parents[i] = Ind(sprand(Bool, windrow, windcol, PB), 0 )
            parents[i].dm[diagind(parents[i].dm)] .=0
            dropzeros!(parents[i].dm)
            evaluate_sw(parents[i], DATASET, off_test, is, js, windrow, windcol)
        end
    end  
    sort!(parents, rev = rv)

    children = Array{Ind}(undef, λ)
    nnz_hist = zeros(epochNum) # number of nonzero elements for best individual in each epoch

    performance = zeros(epochNum, μ) # fitness of parents in each epoch 
    convergence = zeros(epochNum) # how similar are individuals compared to the best individual in the parents  
    greenDot_Num = zeros(Int, 2,epochNum)
    for currentEpoch in 1:epochNum
        if(verbose_)
            @printf "-------current epoch is %d-------\n" currentEpoch
        end 
        for indx in 1:(λ-eli_num)
            if(rand()<0.95-0.95*currentEpoch/epochNum)
                a = LRank(μ, p)
                b = LRank(μ, p)
                while(a==b)
                    b = LRank(μ, p)
                end 

                children[indx] = crossover_deterministic(parents[a], parents[b], windrow, windcol, 1)

            else
                children[indx] = parents[LRank(μ, p)]
            end
        end 

        for indx in 1:(λ-eli_num)
            if(rand() < 0.001+0.999*currentEpoch/epochNum)
                children[indx] = mutate_flip(children[indx])
            end
        end
        children[end-eli_num+1:end] = parents[end-eli_num+1:end]
        for i in 1:λ
            evaluate_sw(children[i], DATASET, off_test, is, js, windrow, windcol)
        end 
        sort!(children, rev = rv)
        parents = children[end-μ+1:end]
        
        
        convergence[currentEpoch] = sum(distance.(Ref(parents[end]), parents, zone_size))
        tmp = map(x->x.error, parents)
        performance[currentEpoch, :] = tmp
        nnz_hist[currentEpoch] = nnz(parents[end].dm)
        

        tmx = green_window.(parents, is, js, windrow, windcol)
        greenDot_Num[1,currentEpoch] = maximum(tmx)
        greenDot_Num[2,currentEpoch] =  green_window(parents[end], is, js, windrow, windcol)
    
        p2 = zeros(μ)
        for i in 1:μ
            p2[i] = green_window(parents[i], is_, js_, window_size_l, window_size_r)
        end 
        maxval = maximum(p2)
        
        pos = [i for (i, x) in enumerate(p2) if x == maxval]
        nnz.(map(x->x.dm, parents))        
        plt1 = plot(1:μ, tmp, label = false, title = "current Epoch "*string(currentEpoch))
        scatter!([pos], [tmp[pos]], label = false)#, ylim = [0, 1])
        plt2 = plot(1:μ, p2, label = false, title = "current Epoch "*string(currentEpoch))
        plot!(1:μ,  nnz.(map(x->x.dm, parents)) , label = false, title = "current Epoch "*string(currentEpoch))
        scatter!([pos], [p2[pos]], label = false, markersize = 3)
        hline!([maxval maxval], label = false)
        hline!([nnz_hist[currentEpoch] nnz_hist[currentEpoch]], label = false)
        hline!([nnz(r_) nnz(r_)], label = false)
        
        display(plot(plt1, plt2, layout=(2,1)))


        if(true)
            @printf "nnz is %d\n" nnz_hist[currentEpoch]  
            @printf "nnz_maxs is %d\n" maxval  
            @printf "best error is: %.4f\n" performance[currentEpoch, end]
        end 

        if(convergence[currentEpoch] > 0.8*μ)
            if(verbose)
                println("!!!RESTART!!!")
            end 
            msk = sample(1:μ, restart_num, replace = false)
            for i in msk
                parents[i] = Ind(sprand(Bool, windrow, windcol, PB),0)
                evaluate_sw(parents[i], DATASET, off_test, is, js, windrow, windcol)
            end 
        end 



    end 

    return parents
end 







p1 = Ind(sprand(Bool, window_size_l, window_size_r, PB),0)
p2 = Ind(sprand(Bool, window_size_l, window_size_r, PB),0)



function roulette_wheel_secltion(population)
    population_fitness = sum([chromosome.error for chromosome in population])
    chromosome_probabilities = [chromosome.error/population_fitness for chromosome in population]
    chromosome_probabilities = 1 .- chromosome_probabilities
    sample(population, Weights(chromosome_probabilities), 2, replace = false)
end 

is_
js_
window_size_l
window_size_r
r_
DATASET_ = r3
r3_mod = DATASET_[1:TIMING, 1:DIM]



