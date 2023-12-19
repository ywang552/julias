tmpData = data_test.-off_test
B = tmpData[1:end-1,:]
A = tmpData[2:end,:]

W = Variable(32,32)
p = minimize(norm(vec(W),2))
p.constraints += A == B*W
p.constraints += diag(W) == 0
solve!(p, SCS.Optimizer; silent_solver = true)


println("true value is: ", rmsd(A, B*reso))
println("true nnz is ", nnz(reso))
x = sparse(W.value)
println("before total error is: ", rmsd(A, B*x), " ", nnz(x))


# x[abs.(x).<=0.1] .=0
# dropzeros!(x)
# nnz(x), rmsd(A,B*x)


# z11 = sparse(B\A)
# # rmsd(B*z,A)
# sum(abs.(z11).<0.01)
# z11[abs.(z11).<0.01] .=0 
# dropzeros!(z11)
# rmsd(B*z11, A)


# # sum(abs.(x).<0.01)
# # x[abs.(x).<0.1] .=0 
# # dropzeros!(x)
# # rmsd(B*x, A)

# # x=z
# th=0.15+0.08
# # th = 0.01
# x_cs = abs.(x).>=th
# x_cg = abs.(x).<th
# dropzeros!(x_cg)
# pb2= (PB*DIM^2-nnz(x_cs))/nnz(x_cg)
# st = rand(DIM,DIM) .< pb2
# function makeParents(d, p, xs, xg)
#     st = rand(d,d).<p
#     x = sparse(xg.* st)
#     x = x+xs
#     # x = xs

#     z = zeros(Bool, d,d)
#     z = sparse(z)
#     z[x.!=0] .= 1
#     xi = Ind(z,0)
#     evaluate_sw(xi, data_test, off_test, is_, js_, window_size_l, window_size_r)
#     return xi
# end 

# ps = Array{Ind}(undef, μ_)
# for i in 1:μ_
#     ps[i] = makeParents(DIM,pb2,x_cs,x_cg)
# end 
# pt="./image/2023.12.4/convex1_th"*string(th)
# g = plotSol(ps[1])
# p = plot(g, title = "Convex child1; th = "*string(th))
# savefig(p,"./image/2023.12.4/convex1_th23")
# g = plotSol(ps[3])
# p = plot(g, title = "Convex child2; th = "*string(th))
# savefig(p,"./image/2023.12.4/convex2_th23")

# x = x_cs
# x.!=0
# z = zeros(32,32)
# z = sparse(z)
# z[x.!=0] .= 1
# xi = Ind(z,0)
# evaluate_sw(xi, data_test, off_test, is_, js_, window_size_l, window_size_r)


# println(nnz(x),", ", nnz(x)-green_window(xi, is_, js_, window_size_l, window_size_r))
# g1 = plotSol(Ind(r_,0))
# g2 = plotSol(xi)

# p = plot(g2, title = "convex using th = "*string(th))
# savefig(p,"./image/2023.12.4/convex_th23")


x=z11
th=0.15+0.08
th = 0.01
x_cs = abs.(x).>=th
x_cg = abs.(x).<th
dropzeros!(x_cg)
pb2= (PB*DIM^2-nnz(x_cs))/nnz(x_cg)
st = rand(DIM,DIM) .< pb2
function makeParents(d, p, xs, xg)
    st = rand(d,d).<p
    x = sparse(xg.* st)
    x = x+xs
    # x = xs

    z = zeros(Bool, d,d)
    z = sparse(z)
    z[x.!=0] .= 1
    xi = Ind(z,0)
    evaluate_sw(xi, data_test, off_test, is_, js_, window_size_l, window_size_r)
    return xi
end 


ps = Array{Ind}(undef, μ_)
for i in 1:μ_
    ps[i] = makeParents(DIM,pb2,x_cs,x_cg)
end 

g = plotSol(ps[1])
p = plot(g, title = "LS child1; th = "*string(th))
savefig(p,"./image/2023.12.4/LS1_th23")
g = plotSol(ps[23])
p = plot(g, title = "LS child2; th = "*string(th))
savefig(p,"./image/2023.12.4/LS2_th23")

x=z11
x.!=0
z = zeros(32,32)
z = sparse(z)
z[x.!=0] .= 1
xi = Ind(z,0)
evaluate_sw(xi, data_test, off_test, is_, js_, window_size_l, window_size_r)


println(nnz(x),", ", nnz(x)-green_window(xi, is_, js_, window_size_l, window_size_r))
g1 = plotSol(Ind(r_,0))
g2 = plotSol(xi)

p = plot(g2, title = "LS using th = "*string(th))
savefig(p,"./image/2023.12.4/LS_th23")







# c = 1
# for i in 0.15:0.01:0.3

#     x_c = abs.(x).>=i

#     sum(abs.(x).<i)

#     # st = rand(DIM,DIM).<0.38
#     # x = sparse(x.* st)
#     # println("removing ", DIM*DIM-sum(nnz(x)), " elements close to 0, leaving ", sum(nnz(x)), " nnz")
#     # println("after total error is: ", rmsd(A, B*x))


#     x_c.!=0
#     z = zeros(32,32)
#     z = sparse(z)
#     z[x_c.!=0] .= 1
#     xi = Ind(z,0)
#     println(nnz(x_c),", ", nnz(x_c)-green_window(xi, is_, js_, window_size_l, window_size_r))
#     a[c] = nnz(x_c)
#     b[c] = nnz(x_c)-green_window(xi, is_, js_, window_size_l, window_size_r)
#     c = c+1
# end 
# # println("number of green dots is ", green_window(xi, is_, js_, window_size_l, window_size_r))
# plot(1:15,a[1:end-1].-a[2:end])
# plot!(1:15,b[1:end-1].-b[2:end])

# # g1 = plotSol(Ind(r_,0))
# # g2 = plotSol(xi)
# # # g3 = plotSol(p[1])
# # plot(g1,g2, layout=(2,1))