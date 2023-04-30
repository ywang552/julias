tt = nnz(r_)
list = findall(!iszero, r_)
list = shuffle(list)
total_trials = 100
h = Array{Float64}(undef, total_trials, 1000)
k = 142

for iter in 1:total_trials
    local c 
    c = 1
    println(iter)
    for indx in k:5:tt
        local z, msk
        z = Ind(spzeros(Bool, window_size_l, window_size_r),0)
        msk = list[1:indx]    
        z.dm[msk] .= 1
        number_to_add = tt-indx
        for indx2 in 1:number_to_add
            row_index = rand(1:window_size_l)
            col_index = rand(1:window_size_r)
            while(z.dm[row_index, col_index] == 1)
                row_index = rand(1:window_size_l)
                col_index = rand(1:window_size_r)
            end 
            z.dm[row_index,col_index] = 1

        end 
        h[iter, c] = evaluate_sw(z, data_test, off_test, is_, js_, window_size_l, window_size_r)
        c = c+1
    end 
end 

# c=size(tt-1500:20:tt)[1]
k_ = Int((tt-k)/5)
h = h[:,1:k_+1]

av = mean(h,dims = 1)'
lw = minimum(h, dims = 1)'
he = maximum(h, dims = 1)'
z = plot(1:k_+1, av)
plot!(1:k_+1, lw)
plot!(1:k_+1, he)
display(z)

plot(1:k_+1, h[4,:])

# plot(1:total_trials, h[:,end])
# hline!([0 ])




# z = Ind(spzeros(Bool, window_size_l, window_size_r),0)
# msk = list[1:375]    
# z.dm[msk] .= 1
# number_to_add = tt-375
# for indx2 in 1:number_to_add
#     row_index = rand(1:window_size_l)
#     col_index = rand(1:window_size_r)
#     while(z.dm[row_index, col_index] == 1)
#         row_index = rand(1:window_size_l)
#         col_index = rand(1:window_size_r)
#     end 
#     z.dm[row_index,col_index] = 1
# end
# evaluate_sw(z, data_test, off_test, is_, js_, window_size_l, window_size_r)

# calculate fitness, the number of connection from 1: total 
# number of connections, the rest is filled with random 
# connections 
# tn = 50
# h = zeros(tn, size(1:tt,1))
# for i in 1:tn
#     println(i)
#     for idx in 1:tt
#         w = idx
#         # println(idx)
#         tc = sample(findall(!iszero, r_), w, replace = false)
#         pt = Ind(spzeros(Bool, window_size_l, window_size_r), 0 )
#         rt = pt.dm
#         rt[tc] .= 1
#         for i in 1:tt-w
#             i, j = rand(1:window_size_l), rand(1:window_size_r)
#             while((i,j) in findnz(pt.dm))
#                 i, j = rand(1:window_size_l), rand(1:window_size_r)
#             end 
#             pt.dm[i,j] = 1
#         end 
#         # plotSol_window(pt, 0, 0, window_size_l, window_size_r)
#         h[i, idx] = evaluate_sw(pt, data_test, off_test, 0, 0, window_size_l, window_size_r)
#     end 
# end 

# av = mean(h,dims = 1)'
# lw = minimum(h, dims = 1)'
# he = maximum(h, dims = 1)'
# plot(1:tt, av)
# plot!(1:tt, lw)
# plot!(1:tt, he)

# pt = Ind(sprand(Bool, window_size_l, window_size_r, PB), 0 )
# xd = evaluate_sw(pt, data_test, off_test, 0, 0, window_size_l, window_size_r)
# hline!([xd, xd])

