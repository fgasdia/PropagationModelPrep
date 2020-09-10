module PropagationModelPrep

using Reexport
using Interpolations

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

include("utils.jl")
include("ComputeJobs.jl")
include("MSISatmosphere.jl")

export ComputeJob, Summit, LocalOMP, Local
export writeshfile, runjob

# Include submodules
include("EMP2D.jl")
include("LWPC.jl")

@reexport using .EMP2D
@reexport using .LWPC

end  # module
