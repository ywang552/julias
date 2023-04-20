## storing resolution 

fp = joinpath(pwd(), "tmp" )
fn = string(runID) * "_reso"
fp_ = (joinpath(fp,fn))
serialize(fp_, r_)

for triNum in 1:100
    local fn, fp_
    a = ES_cheat()
    fn = string(runID)*"_"*string(triNum)*"_solP"
    fp_ = (joinpath(fp,fn))
    serialize(fp_, a)
end 