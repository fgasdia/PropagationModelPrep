module PropagationModelPrep

using Reexport
using Interpolations

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

include("ComputeJobs.jl")
include("MSISatmosphere.jl")

export ComputeJob, Summit
export writeshfile

# Include submodules
include("EMP2D.jl")

@reexport using .EMP2D

end  # module
