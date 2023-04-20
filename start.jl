include("./dependancies.jl")
include("./parameters.jl")
include("./ES.jl")
# include("./solution_analysis.jl")
# include("./result_analysis.jl")

runID = 1 
for run in 1:100
    println("Current RUN is:", runID)
    include("test_generation.jl")
    include("./run.jl")
    runID = runID + 1
end 

