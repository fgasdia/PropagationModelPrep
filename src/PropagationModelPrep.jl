module PropagationModelPrep

using Distributed
using Reexport
using Interpolations
using ProgressMeter

using LongwaveModePropagator
const LMP = LongwaveModePropagator

include("utils.jl")
include("ComputeJobs.jl")

export ComputeJob, Summit, Local, LocalParallel
export writeshfile, runjob

# Include submodules
include("EMP2D.jl")
include("LWPC.jl")

@reexport using .EMP2D
@reexport using .LWPC

end  # module
