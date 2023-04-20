
d = map(x->x.dm, a)
solm = d[:,end]
solutionmap = zeros(window_size_l, window_size_r)
for i in 1:trialNum
    solutionmap = solm[i]+solutionmap
end
solutionmap


b = zeros(trialNum)
for i in 1:trialNum
    b[i] = sum(solutionmap.==i)
end 
b[11]

plot(1:trialNum, b, label = false,  title = "overlapping occurance", xlabel = "hit times", ylabel = "frequency")
hline!([0 0])

msk = (solutionmap.==1)
# # msk = (solutionmap.==2)
# # msk = (solutionmap.==3)
# # msk = (solutionmap.==4)
# # msk = (solutionmap.==5)
# # msk = (solutionmap.==6)
# msk = (solutionmap.==7)
# msk = (solutionmap.==8)

msk = (solutionmap.==6)
filter(x->x[1]!=x[2], findall(!iszero, msk))
wen = zeros(trialNum)
for i in 1:trialNum

    msk = (solutionmap.==i)
    r_t = spzeros(Bool, window_size_l, window_size_r)
    r_t[msk].=1
    dropzeros!(r_t)
    
    
    xd = Ind(r_t,0)
    z1, z2 = calSol_window(xd, is_, js_, window_size_l, window_size_r)
    z1 = z1 == 0 ? 1 : z1
    wen[i] = z2/z1
    # plot(z, title = "hit times = "*string(i))
end 

