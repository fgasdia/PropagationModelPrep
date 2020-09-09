module PropagationModelPrep

using Reexport
using Interpolations

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

include("utils.jl")
include("ComputeJobs.jl")
include("MSISatmosphere.jl")

export ComputeJob, Summit, LocalOMP
export writeshfile, runjob

# Include submodules
include("EMP2D.jl")

@reexport using .EMP2D

end  # module
