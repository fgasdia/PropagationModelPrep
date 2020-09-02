module PropagationModelPrep

using DSP

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

abstract type ComputeJob end

"""
walltime "00:57:00" is 57 minutes
"""
struct Summit <: ComputeJob
    runname::String
    rundir::String
    numnodes::Int
    walltime::String
    exefile::String
end

include("utils.jl")
include("MSISatmosphere.jl")
include("emp2d.jl")

end  # module
