module PropagationModelPrep

using DSP

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

include("utils.jl")
include("MSISatmosphere.jl")

# Include submodules
include("EMP2D.jl")

# using .EMP2D   # not sure if this is needed or not...

end  # module
