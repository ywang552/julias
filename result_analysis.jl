
d = map(x->x.dm, a)
solm = d[:,end]
solutionmap = zeros(window_size_l, window_size_r)
for i in 1:trialNum
    solutionmap = solm[i]+solutionmap
end



b = zeros(trialNum)
for i in 1:trialNum
    b[i] = sum(solutionmap.==i)
end 
b[11]

plot(1:trialNum, b, label = false,  title = "overlapping occurance", xlabel = "hit times", ylabel = "frequency")
hline!([0 0])

msk = (solutionmap.==13)
spp = spzeros(window_size_l, window_size_r)
spp[msk].=1
xd = Ind(spp,0)

plotSol_window(xd, is_, js_, window_size_l, window_size_r)

# # msk = (solutionmap.==2)
# # msk = (solutionmap.==3)
# # msk = (solutionmap.==4)
# # msk = (solutionmap.==5)
# # msk = (solutionmap.==6)
# msk = (solutionmap.==7)
# msk = (solutionmap.==8)


for i in 20:trialNum
    msk = (solutionmap.==i)
    tp = findall(!iszero, msk)
    for j in eachindex(tp)
        k = evaluate_singleGene(tp[j], data_test, off_test)
        println(i,"__", k)
    end 
end 


