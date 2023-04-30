## storing resolution 

# fp = joinpath(pwd(), "tmp" )
# fn = string(runID) * "_reso"
# fp_ = (joinpath(fp,fn))
# serialize(fp_, r_)
r_
EpochNum_p = 150
x = ES_cheat(PB =  0.05)


# x2 = x
# maximum(nnz.([z.dm for z in x]))
# maximum(green_window.(x, is_, js_, window_size_l, window_size_r))
# a = green_window.(x, is_, js_, window_size_l, window_size_r)

# for i in 1:400
#     println(i, "_", green_window(x[i], is_, js_, window_size_l, window_size_r))
# end 

# p1 = zeros(400)
# p2 = zeros(400)
# for i in 1:400
#     p1[i] = evaluate_sw(x[i], data_test, off_test, is_,js_,window_size_l, window_size_r)
#     p2[i] = green_window(x[i], is_, js_, window_size_l, window_size_r)
# end 
# maxval = maximum(p2)

# pos = [i for (i, x) in enumerate(p2) if x == maxval]

# plt1 = plot(1:400, p1, label = false)
# scatter!([pos], [p1[pos]], label = false)
# plt2 = plot(1:400, p2, label = false)
# scatter!([pos], [p2[pos]], label = false, markersize = 3)
# plot(plt1, plt2, layout=(2,1))













# savefig("./image/issue with incremental")

# sd = Array{Ind}(undef, μ_)
# for i in 0.08:0.05:0.38
#     global sd
#     if(i == 0.08)
#         sd = ES_cheat(PB =  i)
#     else
#         seeded_ = true
#         seed_ = sd
#         sd = ES_cheat(seeded = seeded_, seed = seed_, PB = i)
#     end 
# end 

# for triNum in 1:100
#     local fn, fp_
#     a = ES_cheat()
#     fn = string(runID)*"_"*string(triNum)*"_solP"
#     fp_ = (joinpath(fp,fn))
#     serialize(fp_, a)
# end 