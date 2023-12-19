tmpData = data_test.-off_test
# B = r1[1:end-1,:]
# A = r1[2:end,:]
B = tmpData[1:end-1,:]
A = tmpData[2:end,:]
x = sparse(B\A)
x[abs.(x).<0.15].=0
dropzeros!(x)
x
r_

function countSol(ptr)
    m = ptr.dm
    c1 = 0
    c2 = 0 
    c3 = 0
    for i in eachindex(r_)
        if(m[i] == 1 && r_[i] == 1)
            c1 = c1+1
        end 
        if(m[i] == 1 && r_[i] == 0)

            c2 = c2+1
        end 
        if(m[i] == 0 && r_[i] == 1 )
    
            c3 = c3+1
        end 
    end 
    return c1,c2,c3
end 
z = zeros(DIM,DIM)
z = sparse(z)
z[x.!=0] .= 1
xi = Ind(z,0)
plotSol(Ind(r_,0))
plotSol(xi)
# savefig("xd")
# countSol(Ind(r_,0))
# countSol(xi)

# r_
# C = B*x
# rmsd(B*x,A)
# g= 3
# plot(A[:,g])
# plot!(C[:,g])

# z = zeros(DIM,DIM)
# z = sparse(z)
# z[x.!=0] .= 1
# dropzeros!(z)
# p1 = Ind(z,0)
# plotSol(Ind(r_,0))
# plotSol(p1)

# B = r1[1:end-1,:]
# A = r1[2:end,:]
plotSol(xi)
m = xi.dm
for j in 1:DIM
    pos = j
    notTerminated = true
    idx = findall(x->x==1, m[:,pos])
    x_ = B[:,idx]\A[:,pos]
    c2 = rmsd(B[:,idx]*x_,A[:,pos]) .>0.0005

    while(notTerminated && c2)
        x_ = B[:,idx]\A[:,pos]
        rtA = ones(DIM)
        rtmin = rmsd(B[:,idx]*x_,A[:,pos])

        for i in 1:DIM
            push!(idx,i)
            x_ = B[:,idx]\A[:,pos]
            rt = rmsd(B[:,idx]*x_,A[:,pos])
            if(rt < rtmin)
                rtA[i] = rt
            end 
            idx = idx[1:end-1]
        end
        mn = findmin(rtA)
        push!(idx,mn[2])
        if(abs(mn[1]) < 0.0001)
            notTerminated = false
        end 
    end 
    println("---------", j,"-----------")
    println(sort(idx))
    println(sort(findall(x->x==1,r_[:,pos])))
end 
