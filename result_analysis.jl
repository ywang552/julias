fp = joinpath(pwd(), "tmp")
hitTime = zeros(10, 100)
accuracy = zeros(10, 100)
for runID in 1:10 
    global fp
    local r_, fn
    fn = string(runID)*"_reso"
    r_ = deserialize(joinpath(fp, fn))
    p_best = Array{Ind}(undef, 100)
    for i in 1:100
        local fn
        fn = string(runID)*"_"*string(i)*"_solP"
        parents = deserialize(joinpath(fp,fn))
        p_best[i] = parents[end]
    end 
    dms = map(x->x.dm, p_best)
    tp1 = spzeros(window_size_l, window_size_r)
    for v in eachindex(dms)
        local msk
        msk = diagind(dms[v])
        dms[v][msk].=0
        dropzeros!(dms[v])
    end 


    for i in 1:100
        tp1 = tp1 + dms[i]
    end 

    
    for i  in 1:100
        global  hitTime
        local msk 
        tmpInd = Ind(spzeros(Bool, window_size_l, window_size_r), 0)
        msk = tp1 .== i 
        hitTime[runID, i] = sum(msk)
        tmpInd.dm[msk] .=1 
        c1,c2 = calSol_window(tmpInd, is_, js_, window_size_l, window_size_r)
        c2d = maximum([1, c2])
        c1d = maximum([1, c1])
        accuracy[runID, i] = c1
    end 

end 

z = 1
pl1 = plot(1:100, hitTime[z,:], legend = false )
plot!(1:100, accuracy[z,:], legend = false)






