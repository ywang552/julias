test_dim = DIM
vc = 37
timet_test = TIMING
off_test = rand(Uniform(0, 10), 1,test_dim, )
reso = sprand(test_dim, test_dim,PB, n->rand(Uniform(-0.5,0.5),n)) 
reso[diagind(reso)] .= 0
r_ = spzeros(Bool, DIM,DIM)
r_[reso.!=0] .= 1

init_test = rand(1, test_dim) 
init_test
data_test = zeros(timet_test, test_dim)
data_test[1,:] = init_test
ubt = off_test.+0.5
lbt = off_test.-0.5

for i in 2:timet_test
    data_test[i, :] = data_test[i-1,:]'*reso 
end

data_test = data_test.+off_test
pl = plot(1:vc, data_test[1:vc, 1:10], legend = false)
display(pl)
sleep(2)
@printf "reso size: %d\n" nnz(reso)

