
fp = joinpath(pwd(), "data", "259584", "3_hist")
a = deserialize(fp)


histbook = a[1]
nnzbook = a[2]

u = plot(1:5000, histbook[1:5000, 1], xaxis=:log)
plot!(1:5000, histbook[1:5000,8],xaxis=:log)
b = plot(1:5000, nnzbook[1:5000],xaxis=:log)
display(plot(u,b,layout=(2,1)))

p_t = Array{Ind}(undef, 10, 16)
for i in 1:10 
    local fp_
    fp_ = joinpath(pwd(), "data", "259584", string(i)*"_models")
    p_t[i,:] = deserialize(fp_)
end 

r = spzeros(128,128)

for i in 1:10
    global r
    r = r + p_t[i,1].dm
end 
s
println(sum(r.==1))
println(sum(r.==2))
println(sum(r.==3))

